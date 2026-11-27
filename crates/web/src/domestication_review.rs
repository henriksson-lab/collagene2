//! Exact-event installation of independently checked substitution intent.
//! Opening a candidate does not retain positive reaction verification.
use collagene_core::{
    Molecule, domestication as design, domestication_preparation as physical,
    history::DocumentHistory,
    run_lineage::{EventBinding, InputEventBinding},
    transfer::Candidate,
    workspace::{Workspace, WorkspaceAttachment},
};
#[derive(Clone)]
pub struct Open {
    pub source: EventBinding,
    pub origin: Option<crate::run_origin_model::PinnedOrigin>,
    pub epoch: u64,
    pub request: design::Request,
    pub selection: design::Selection,
    pub check: design::CheckReport,
    pub target: Molecule,
    pub prepared: Option<physical::Prepared>,
}
pub struct Staged {
    pub workspace: Workspace,
    pub target: EventBinding,
    pub inputs: Vec<InputEventBinding>,
    pub candidate: Option<Candidate>,
}
fn fingerprint(value: &impl serde::Serialize) -> Result<String, collagene_core::Error> {
    struct Writer {
        hash: blake3::Hasher,
        bytes: usize,
    }
    impl std::io::Write for Writer {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.bytes = self
                .bytes
                .checked_add(bytes.len())
                .filter(|n| *n <= 32 * 1024 * 1024)
                .ok_or_else(|| std::io::Error::other("domestication comparison bound"))?;
            self.hash.update(bytes);
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut writer = Writer {
        hash: blake3::Hasher::new(),
        bytes: 0,
    };
    serde_json::to_writer(&mut writer, value).map_err(|_| collagene_core::Error::TooLarge)?;
    Ok(writer.hash.finalize().to_hex().to_string())
}
/// Attest an independently captured material selection against the exact native
/// event and complete source snapshot. No replay or material restoration occurs.
pub fn validate_origin(
    frame: &Workspace,
    event: &EventBinding,
    source: &Molecule,
    origin: Option<&crate::run_origin_model::PinnedOrigin>,
) -> Result<(), String> {
    if let Some(origin) = origin {
        if origin.event() != event {
            return Err("Saved material belongs to another source event.".into());
        }
        origin
            .as_ref()
            .validate_for_workspace(frame, event)
            .map_err(|e| e.to_string())?;
        if frame
            .revisions
            .get(&event.revision)
            .ok_or("Native source revision is absent")?
            != source
        {
            return Err("Saved material source snapshot changed.".into());
        }
    }
    Ok(())
}
pub fn stage(frame: &Workspace, opened: &Open) -> Result<Staged, String> {
    frame.validate().map_err(|e| e.to_string())?;
    let docs = DocumentHistory::restore_all(frame).map_err(|e| e.to_string())?;
    let original = docs
        .iter()
        .find(|d| d.id == opened.source.document_id)
        .ok_or("The original domestication source is unavailable.")?;
    let source = original
        .revision_molecule(&opened.source.event_id)
        .map_err(|e| e.to_string())?;
    if source.revision() != opened.source.revision
        || opened.request.source_revision != opened.source.revision
    {
        return Err("Domestication must retain its exact original source event.".into());
    }
    validate_origin(frame, &opened.source, &source, opened.origin.as_ref())?;
    let checked =
        design::check(&source, &opened.request, &opened.selection).map_err(|e| e.to_string())?;
    opened
        .check
        .validate_binding(&checked.report().context_revision)
        .map_err(|e| e.to_string())?;
    if opened.target != *checked.target()
        || fingerprint(&opened.check).map_err(|e| e.to_string())?
            != fingerprint(checked.report()).map_err(|e| e.to_string())?
    {
        return Err("The intended substitutions, features or checked target changed.".into());
    }
    let candidate = if let Some(prepared) = &opened.prepared {
        let input = collagene_core::verification::PlanInput {
            id: physical::INPUT_ID.into(),
            revision: source.revision(),
            molecule: source.clone(),
        };
        if fingerprint(&prepared.request.origin).map_err(|e| e.to_string())?
            != fingerprint(&opened.origin.as_ref().map(|origin| origin.as_ref()))
                .map_err(|e| e.to_string())?
        {
            return Err(
                "The preparation material origin differs from the independently selected source."
                    .into(),
            );
        }
        let expected = physical::revision(&input, &prepared.request, Some(frame))
            .map_err(|e| e.to_string())?;
        prepared
            .validate_binding(&expected)
            .map_err(|e| e.to_string())?;
        if fingerprint(&prepared.request.design).map_err(|e| e.to_string())?
            != fingerprint(&opened.request).map_err(|e| e.to_string())?
            || fingerprint(&prepared.request.selection).map_err(|e| e.to_string())?
                != fingerprint(&opened.selection).map_err(|e| e.to_string())?
            || prepared.target != opened.target
            || prepared.plan.inputs.len() != 1
            || prepared.plan.inputs[0].id != physical::INPUT_ID
            || prepared.plan.inputs[0].molecule != source
        {
            return Err("The proposed physical preparation differs from the checked original source and independent intent.".into());
        }
        Some(Candidate::new(prepared.plan.clone()))
    } else {
        None
    };
    let mut staged = frame.clone();
    let existing = docs.iter().find(|d| {
        d.head_revision() == opened.target.revision()
            && d.summaries()
                .iter()
                .any(|event| event.current && event.parents == [opened.source.event_id.clone()])
    });
    let target = if let Some(doc) = existing {
        EventBinding {
            document_id: doc.id.clone(),
            event_id: doc.head_event().into(),
            revision: doc.head_revision().into(),
        }
    } else {
        if docs.len() >= collagene_core::workspace::MAX_DOCUMENTS {
            return Err("The workspace document limit has been reached.".into());
        }
        let prefix = format!("domestication-{}", &opened.check.context_revision[..16]);
        let id = (1..=1001)
            .map(|n| format!("{prefix}-{n}"))
            .find(|id| !docs.iter().any(|d| d.id == *id))
            .ok_or("No intended-reference identity is available.")?;
        let derived=original.derive_event(&opened.source.event_id,id,opened.target.clone(),"Checked source-coordinate substitution intent; reference only, no experimental outcome established").map_err(|e|e.to_string())?;
        let binding = EventBinding {
            document_id: derived.id.clone(),
            event_id: derived.head_event().into(),
            revision: derived.head_revision().into(),
        };
        derived.contribute(&mut staged).map_err(|e| e.to_string())?;
        binding
    };
    #[derive(serde::Serialize)]
    struct ReviewEvidence<'a> {
        version: u16,
        source: &'a EventBinding,
        origin: Option<&'a collagene_core::run_material::MaterialOrigin>,
        target_event: &'a EventBinding,
        check: &'a design::CheckReport,
        // Always serialized, including explicit null for intent-only staging.
        prepared: Option<&'a physical::Prepared>,
        authority: &'static str,
    }
    struct EvidenceBytes(Vec<u8>);
    impl std::io::Write for EvidenceBytes {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            if bytes.len() > collagene_core::workspace::MAX_ENTRY_BYTES.saturating_sub(self.0.len())
            {
                return Err(std::io::Error::other(
                    "Domestication review exceeds its native attachment bound",
                ));
            }
            self.0.extend_from_slice(bytes);
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut output = EvidenceBytes(Vec::new());
    serde_json::to_writer(&mut output,&ReviewEvidence{
        version:physical::REVIEW_EVIDENCE_VERSION,
        source:&opened.source,
        origin:opened.origin.as_ref().map(|origin|origin.as_ref()),
        target_event:&target,
        check:checked.report(),
        prepared:opened.prepared.as_ref(),
        authority:"Historical documentary source/intent/preparation snapshot. Recorded primer readiness and verification are not restored as current authority; fresh reconstruction is required.",
    }).map_err(|e|e.to_string())?;
    let evidence = output.0;
    let hash = blake3::hash(&evidence).to_hex().to_string();
    staged
        .attachments
        .entry(hash)
        .or_insert(WorkspaceAttachment {
            kind: physical::REVIEW_EVIDENCE_KIND.into(),
            name: "domestication-intent.json".into(),
            media_type: "application/json".into(),
            bytes: evidence,
        });
    staged.active_document = target.document_id.clone();
    let inputs = vec![InputEventBinding {
        input_id: physical::INPUT_ID.into(),
        event: opened.source.clone(),
    }];
    if let Some(candidate) = &candidate {
        staged =
            crate::assembly_review::stage_existing(&staged, candidate, &target, &inputs)?.workspace;
    } else {
        staged.active_candidate = None;
    }
    staged.validate().map_err(|e| e.to_string())?;
    Ok(Staged {
        workspace: staged,
        target,
        inputs,
        candidate,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::{Feature, Location, Topology, catalog_digest, codon_design};
    fn fixture() -> (Workspace, Open) {
        let mut source = Molecule::dna(
            "original annotated plasmid",
            "ATGGGTCTCTAACCAGTATCGACGATATCACGTAACCTGCACTAGTTCAGGACTGACATCGAACGTTAGTC",
            Topology::Circular,
        )
        .unwrap();
        source.features = vec![
            Feature {
                kind: "CDS".into(),
                location: Location::range(0, 12),
                qualifiers: vec![("gene".into(), Some("unchanged protein".into()))],
                needs_review: false,
            },
            Feature {
                kind: "misc_feature".into(),
                location: Location::range(6, 7),
                qualifiers: vec![(
                    "note".into(),
                    Some("one unchanged base between substitutions".into()),
                )],
                needs_review: false,
            },
        ];
        let frame =
            Workspace::from_documents("substitution ancestry", vec![source.clone()], 0).unwrap();
        let request = design::Request {
            source_revision: source.revision(),
            selection: catalog_digest::select(
                &catalog_digest::builtin_catalog().unwrap(),
                &["bsai".into()],
            )
            .unwrap(),
            scope: vec![design::SourceRange {
                start: 0,
                end: source.sequence.len(),
            }],
            retained_sites: vec![],
            coding_features: vec![design::feature_pin(&source, 0).unwrap()],
            protected: vec![],
            options: codon_design::DesignOptions::default(),
            limits: design::Limits::default(),
        };
        let patches = vec![
            design::Patch {
                position: 5,
                expected: "T".into(),
                replacement: "A".into(),
            },
            design::Patch {
                position: 8,
                expected: "C".into(),
                replacement: "T".into(),
            },
        ];
        let selection = design::Selection {
            context_revision: design::propose_revision(&source, &request).unwrap(),
            candidate_id: design::candidate_id(&request, &patches).unwrap(),
            patches,
        };
        let checked = design::check(&source, &request, &selection).unwrap();
        let original = &frame.documents[0];
        let source = EventBinding {
            document_id: original.id.clone(),
            event_id: original.head_event.clone().unwrap(),
            revision: original.head.clone(),
        };
        let opened = Open {
            origin: None,
            source,
            epoch: 4,
            request,
            selection,
            check: checked.report().clone(),
            target: checked.target().clone(),
            prepared: None,
        };
        (frame, opened)
    }
    #[test]
    fn reference_open_preserves_intervening_feature_coordinates_and_all_source_heads() {
        let (frame, opened) = fixture();
        let before = serde_json::to_vec(&frame).unwrap();
        let original = frame.head(&frame.documents[0].id).unwrap();
        let staged = stage(&frame, &opened).unwrap();
        let target = staged.workspace.head(&staged.target.document_id).unwrap();
        assert_eq!(&target.sequence.as_str()[..12], "ATGGGACTTTAA");
        assert_eq!(target.features[1], original.features[1]);
        assert_eq!(target.features[0].location, original.features[0].location);
        assert_eq!(
            target.features[0].qualifiers,
            original.features[0].qualifiers
        );
        assert!(target.features[0].needs_review);
        assert_eq!(
            staged.workspace.documents[0].head_event,
            frame.documents[0].head_event
        );
        let evidence = staged
            .workspace
            .attachments
            .values()
            .find(|a| a.kind == physical::REVIEW_EVIDENCE_KIND)
            .unwrap();
        let record: serde_json::Value = serde_json::from_slice(&evidence.bytes).unwrap();
        assert_eq!(record["version"], physical::REVIEW_EVIDENCE_VERSION);
        assert!(
            record
                .get("prepared")
                .expect("required nullable preparation")
                .is_null()
        );
        assert_eq!(
            record["source"],
            serde_json::to_value(&opened.source).unwrap()
        );
        assert_eq!(
            record["target_event"],
            serde_json::to_value(&staged.target).unwrap()
        );
        assert!(
            record
                .get("origin")
                .expect("required nullable origin")
                .is_null()
        );
        assert!(staged.candidate.is_none());
        assert!(staged.workspace.active_candidate.is_none());
        let again = stage(&staged.workspace, &opened).unwrap();
        assert_eq!(again.target, staged.target);
        assert_eq!(again.workspace.documents.len(), 2);
        assert_eq!(serde_json::to_vec(&frame).unwrap(), before);
    }
    #[test]
    fn wrong_event_forged_target_or_feature_qualifier_cannot_install() {
        let (frame, opened) = fixture();
        let mut forged = opened.clone();
        forged.target.features[1].location = Location::range(5, 9);
        assert!(stage(&frame, &forged).is_err());
        let mut forged = opened.clone();
        forged.target.features[0].qualifiers.clear();
        assert!(stage(&frame, &forged).is_err());
        let mut wrong = opened.clone();
        wrong.source.event_id = "event-missing".into();
        assert!(stage(&frame, &wrong).is_err());
        let mut changed = opened;
        changed.selection.patches[0].replacement = "C".into();
        assert!(stage(&frame, &changed).is_err());
    }
    #[test]
    fn physical_plan_open_keeps_original_event_and_independent_target_without_restoring_proof() {
        let (frame, mut opened) = fixture();
        let source = frame.head(&frame.documents[0].id).unwrap().clone();
        let input = collagene_core::verification::PlanInput {
            id: physical::INPUT_ID.into(),
            revision: source.revision(),
            molecule: source,
        };
        opened.prepared = Some(
            physical::materialize(
                &input,
                &physical::Request {
                    origin: None,
                    design: opened.request.clone(),
                    selection: opened.selection.clone(),
                    route: physical::Route::InversePcr {
                        edit_start: 5,
                        edit_end: 9,
                        forward_length: 8,
                        reverse_length: 8,
                        reverse_tail_bases: 0,
                    },
                },
                None,
            )
            .unwrap(),
        );
        let staged = stage(&frame, &opened).unwrap();
        let candidate = staged.candidate.unwrap();
        assert_eq!(candidate.plan.target_revision, opened.target.revision());
        assert_eq!(candidate.plan.inputs[0].molecule, input.molecule);
        assert_eq!(staged.inputs[0].event, opened.source);
        assert_eq!(
            staged.workspace.active_candidate,
            Some(candidate.plan_revision)
        );
        assert!(
            staged
                .workspace
                .attachments
                .values()
                .all(|a| a.kind != "cloning-run-v1")
        );
    }
    #[test]
    fn equal_dna_aliases_keep_distinct_source_ancestry() {
        let (mut frame, opened) = fixture();
        let docs = DocumentHistory::restore_all(&frame).unwrap();
        let alias = docs[0]
            .fork_event(&opened.source.event_id, "equal-dna-other-source")
            .unwrap();
        let mut other = opened.clone();
        other.source = EventBinding {
            document_id: alias.id.clone(),
            event_id: alias.head_event().into(),
            revision: alias.head_revision().into(),
        };
        alias.contribute(&mut frame).unwrap();
        frame.validate().unwrap();
        let first = stage(&frame, &opened).unwrap();
        let second = stage(&first.workspace, &other).unwrap();
        assert_ne!(first.target.document_id, second.target.document_id);
        let event = second
            .workspace
            .history
            .iter()
            .find(|e| e.id == second.target.event_id)
            .unwrap();
        assert_eq!(event.parents, [other.source.event_id]);
    }
    #[test]
    fn overlap_extension_linear_and_circle_stage_real_canonical_products_and_archive_sources() {
        use collagene_core::{
            cloning_outputs, codon_design,
            verification::{self, Status},
        };
        for topology in [Topology::Linear, Topology::Circular] {
            let source = Molecule::dna(
                "Original OE parent",
                "AGTCCGATGCTACGTAGAATTCTTCAGACCTGTAGATCGTCA",
                topology,
            )
            .unwrap();
            let frame =
                Workspace::from_documents("OE domestication", vec![source.clone()], 0).unwrap();
            let doc = &frame.documents[0];
            let event = EventBinding {
                document_id: doc.id.clone(),
                event_id: doc.head_event.clone().unwrap(),
                revision: doc.head.clone(),
            };
            let request = design::Request {
                source_revision: source.revision(),
                selection: collagene_core::catalog_digest::select(
                    &collagene_core::catalog_digest::builtin_catalog().unwrap(),
                    &["ecori".into()],
                )
                .unwrap(),
                scope: vec![design::SourceRange {
                    start: 0,
                    end: source.sequence.len(),
                }],
                retained_sites: vec![],
                coding_features: vec![],
                protected: vec![],
                options: codon_design::DesignOptions::default(),
                limits: design::Limits::default(),
            };
            let patches = vec![design::Patch {
                position: 18,
                expected: "A".into(),
                replacement: "C".into(),
            }];
            let selection = design::Selection {
                context_revision: design::propose_revision(&source, &request).unwrap(),
                candidate_id: design::candidate_id(&request, &patches).unwrap(),
                patches,
            };
            let checked = design::check(&source, &request, &selection).unwrap();
            assert_eq!(
                checked.target().sequence.as_str(),
                "AGTCCGATGCTACGTAGACTTCTTCAGACCTGTAGATCGTCA"
            );
            let mut options = collagene_core::overlap_extension_primer_design::Options::default();
            options.primer_options.min_length = 8;
            options.primer_options.max_length = 8;
            options.primer_options.min_gc_fraction = 0.;
            options.primer_options.max_gc_fraction = 1.;
            options.primer_options.max_tm_difference_c = 100.;
            options.primer_options.max_homopolymer = 200;
            options.primer_options.max_three_prime_complementarity = 200;
            options.homology_flank = 8;
            options.min_paired_bases = 8;
            let prepared = physical::materialize(
                &collagene_core::verification::PlanInput {
                    id: physical::INPUT_ID.into(),
                    revision: source.revision(),
                    molecule: source.clone(),
                },
                &physical::Request {
                    origin: None,
                    design: request.clone(),
                    selection: selection.clone(),
                    route: physical::Route::OverlapExtension { options },
                },
                None,
            )
            .unwrap();
            assert!(prepared.proposal_ready, "{:?}", prepared.proposal_findings);
            assert_eq!(
                prepared.plan.steps.len(),
                if topology == Topology::Circular { 6 } else { 4 }
            );
            assert_eq!(
                prepared.verification.status,
                Status::ConsistentUnderAssumptions,
                "{:?}",
                prepared.verification.findings
            );
            let opened = Open {
                origin: None,
                source: event.clone(),
                epoch: 9,
                request: request.clone(),
                selection: selection.clone(),
                check: checked.report().clone(),
                target: checked.target().clone(),
                prepared: Some(prepared),
            };
            let staged = stage(&frame, &opened).unwrap();
            let candidate = staged.candidate.unwrap();
            assert_eq!(staged.workspace.documents[0].head, frame.documents[0].head);
            assert_eq!(
                staged.workspace.documents[0].head_event,
                frame.documents[0].head_event
            );
            assert_eq!(
                verification::verify(&candidate.plan, &opened.target).status,
                Status::ConsistentUnderAssumptions
            );
            let staged_target = staged.target.clone();
            let committed = cloning_outputs::finalize(
                &staged.workspace,
                &candidate.plan,
                staged.target,
                staged.inputs,
                &cloning_outputs::Options {
                    name: "Reviewed OE product".into(),
                    show_prepared: true,
                    save_primers: true,
                    allow_incomplete: false,
                },
            )
            .unwrap();
            assert!(committed.receipt.primer_attachment.is_some());
            let reopened = collagene_core::workspace::load_workspace(
                &collagene_core::workspace::save_workspace(&committed.workspace).unwrap(),
            )
            .unwrap();
            let evidence = reopened
                .attachments
                .values()
                .find(|a| a.kind == physical::REVIEW_EVIDENCE_KIND)
                .unwrap();
            let record: serde_json::Value = serde_json::from_slice(&evidence.bytes).unwrap();
            assert_eq!(record["version"], physical::REVIEW_EVIDENCE_VERSION);
            assert_eq!(record["source"], serde_json::to_value(&event).unwrap());
            assert_eq!(
                record["target_event"],
                serde_json::to_value(&staged_target).unwrap()
            );
            let historical: physical::Prepared = serde_json::from_value(
                record
                    .get("prepared")
                    .expect("required preparation field")
                    .clone(),
            )
            .unwrap();
            let original = opened.prepared.as_ref().unwrap();
            historical
                .validate_binding(&original.context_revision)
                .unwrap();
            assert_eq!(
                serde_json::to_vec(&historical).unwrap(),
                serde_json::to_vec(original).unwrap()
            );
            assert_eq!(historical.proposal_ready, original.proposal_ready);
            assert_eq!(historical.proposal_findings, original.proposal_findings);
            let final_product = &reopened.revisions[&committed.receipt.construct.revision];
            assert_eq!(final_product.sequence, opened.target.sequence);
            assert_eq!(final_product.topology, topology);
            assert_eq!(
                reopened
                    .documents
                    .iter()
                    .find(|d| d.id == event.document_id)
                    .unwrap()
                    .head_event
                    .as_ref(),
                Some(&event.event_id)
            );
            let mut stale = opened.clone();
            stale.source.event_id = "unretained-equal-DNA-event".into();
            assert!(stage(&frame, &stale).is_err());
            let mut protected = request.clone();
            protected
                .protected
                .push(design::SourceRange { start: 18, end: 19 });
            let mut protected_selection = selection.clone();
            protected_selection.context_revision =
                design::propose_revision(&source, &protected).unwrap();
            protected_selection.candidate_id =
                design::candidate_id(&protected, &protected_selection.patches).unwrap();
            assert!(design::check(&source, &protected, &protected_selection).is_err());
            let mut wrong = opened;
            wrong.target.sequence = collagene_core::Sequence::dna("TTTTCCCCAAAAGGGG").unwrap();
            assert!(stage(&frame, &wrong).is_err());
        }
    }
    #[test]
    fn diagnostic_nonready_preparation_remains_documentary_after_native_restart() {
        let source = Molecule::dna(
            "Repeated original endpoint",
            "AGTCCGATAGTCCGATGCTACGTAGAATTCCAGACCTGTAGATCGTCA",
            Topology::Linear,
        )
        .unwrap();
        let frame =
            Workspace::from_documents("Diagnostic OE history", vec![source.clone()], 0).unwrap();
        let mut form = crate::domestication_form::Form::new(&source, None);
        form.enzymes = "ecori".into();
        let request = form.request(&source).unwrap();
        let patches = vec![design::Patch {
            position: 26,
            expected: "A".into(),
            replacement: "C".into(),
        }];
        let selection = design::Selection {
            context_revision: design::propose_revision(&source, &request).unwrap(),
            candidate_id: design::candidate_id(&request, &patches).unwrap(),
            patches,
        };
        let checked = design::check(&source, &request, &selection).unwrap();
        let mut options = collagene_core::overlap_extension_primer_design::Options::default();
        options.primer_options.min_length = 8;
        options.primer_options.max_length = 8;
        options.primer_options.min_gc_fraction = 0.;
        options.primer_options.max_gc_fraction = 1.;
        options.primer_options.max_tm_difference_c = 100.;
        options.primer_options.max_homopolymer = 200;
        options.primer_options.max_three_prime_complementarity = 200;
        options.homology_flank = 8;
        options.min_paired_bases = 8;
        let prepared = physical::materialize(
            &collagene_core::verification::PlanInput {
                id: physical::INPUT_ID.into(),
                revision: source.revision(),
                molecule: source,
            },
            &physical::Request {
                origin: None,
                design: request.clone(),
                selection: selection.clone(),
                route: physical::Route::OverlapExtension { options },
            },
            None,
        )
        .unwrap();
        assert!(!prepared.proposal_ready);
        assert!(!prepared.proposal_findings.is_empty());
        let doc = &frame.documents[0];
        let event = EventBinding {
            document_id: doc.id.clone(),
            event_id: doc.head_event.clone().unwrap(),
            revision: doc.head.clone(),
        };
        let opened = Open {
            origin: None,
            source: event,
            epoch: 3,
            request,
            selection,
            check: checked.report().clone(),
            target: checked.target().clone(),
            prepared: Some(prepared),
        };
        let staged = stage(&frame, &opened).unwrap();
        assert!(staged.candidate.is_some());
        let reopened = collagene_core::workspace::load_workspace(
            &collagene_core::workspace::save_workspace(&staged.workspace).unwrap(),
        )
        .unwrap();
        let asset = reopened
            .attachments
            .values()
            .find(|a| a.kind == physical::REVIEW_EVIDENCE_KIND)
            .unwrap();
        let record: serde_json::Value = serde_json::from_slice(&asset.bytes).unwrap();
        let historical: physical::Prepared =
            serde_json::from_value(record.get("prepared").unwrap().clone()).unwrap();
        let original = opened.prepared.unwrap();
        historical
            .validate_binding(&original.context_revision)
            .unwrap();
        assert!(!historical.proposal_ready);
        assert_eq!(historical.proposal_findings, original.proposal_findings);
        assert_eq!(historical.verification.status, original.verification.status);
        assert!(historical.verification.products.is_empty());
        assert_eq!(
            reopened.documents[0].head_event,
            frame.documents[0].head_event
        );
    }
    #[test]
    fn derived_source_open_retains_exact_origin_and_only_new_oligos_after_restart() {
        use crate::run_origin_model::PinnedOrigin;
        use collagene_core::{
            analysis::Primer,
            chemistry::Strand,
            material_derivation as native,
            run_material::MaterialOrigin,
            verification::{self, PlanInput},
            workspace,
        };
        let source = Molecule::dna(
            "Seed reference",
            "AGTCCGATGCTACGTAGAATTCCAGACCTGTAGATCGTCA",
            Topology::Linear,
        )
        .unwrap();
        let frame =
            Workspace::from_documents("Saved domestication source", vec![source.clone()], 0)
                .unwrap();
        let doc = &frame.documents[0];
        let seed = EventBinding {
            document_id: doc.id.clone(),
            event_id: doc.head_event.clone().unwrap(),
            revision: doc.head.clone(),
        };
        let forward =
            Primer::new("Seed stock forward", &source.sequence.as_str()[..8], "").unwrap();
        let reverse = Primer::new(
            "Seed stock reverse",
            collagene_core::Sequence::dna(&source.sequence.as_str()[source.sequence.len() - 8..])
                .unwrap()
                .reverse_complement()
                .unwrap()
                .as_str(),
            "",
        )
        .unwrap();
        let source_origin = native::SourceOrigin::OriginalReference;
        let previews =
            native::preview_pcr(&source, &source_origin, &forward, &reverse, Strand::Top).unwrap();
        assert_eq!(previews.len(), 1);
        let installed = native::prepare_install(
            &frame,
            &native::Request {
                source: seed,
                origin: source_origin.clone(),
                name: "Saved PCR stock".into(),
                save_primers: true,
                derivation: native::Derivation::ExactPcr {
                    template_strand: Strand::Top,
                    selection: native::PcrSelection {
                        selected_index: 0,
                        expected_preview_revision: previews[0].molecule.revision(),
                        geometry: previews[0].geometry.clone(),
                    },
                    forward,
                    reverse,
                },
            },
        )
        .unwrap();
        let frame =
            workspace::load_workspace(&workspace::save_workspace(&installed.workspace).unwrap())
                .unwrap();
        let event = installed.receipt.product;
        let origin = MaterialOrigin::Derived(native::capsule(&frame, &event).unwrap());
        let pin = PinnedOrigin::from_material(origin.clone()).unwrap();
        let source = frame.revisions[&event.revision].clone();
        let mut form = crate::domestication_form::Form::new(&source, None);
        form.enzymes = "ecori".into();
        let request = form.request(&source).unwrap();
        let patches = vec![design::Patch {
            position: 18,
            expected: "A".into(),
            replacement: "C".into(),
        }];
        let selection = design::Selection {
            context_revision: design::propose_revision(&source, &request).unwrap(),
            candidate_id: design::candidate_id(&request, &patches).unwrap(),
            patches,
        };
        let checked = design::check(&source, &request, &selection).unwrap();
        assert_eq!(
            checked.target().sequence.as_str(),
            "AGTCCGATGCTACGTAGACTTCCAGACCTGTAGATCGTCA"
        );
        let mut options = collagene_core::overlap_extension_primer_design::Options::default();
        options.primer_options.min_length = 8;
        options.primer_options.max_length = 8;
        options.primer_options.min_gc_fraction = 0.;
        options.primer_options.max_gc_fraction = 1.;
        options.primer_options.max_tm_difference_c = 100.;
        options.primer_options.max_homopolymer = 200;
        options.primer_options.max_three_prime_complementarity = 200;
        options.homology_flank = 8;
        options.min_paired_bases = 8;
        let input = PlanInput {
            id: physical::INPUT_ID.into(),
            revision: source.revision(),
            molecule: source.clone(),
        };
        let mut native_request = physical::Request {
            design: request.clone(),
            selection: selection.clone(),
            route: physical::Route::OverlapExtension { options },
            origin: Some(origin.clone()),
        };
        let prepared = physical::materialize(&input, &native_request, Some(&frame)).unwrap();
        assert!(prepared.proposal_ready, "{:?}", prepared.proposal_findings);
        let oligos = collagene_core::plan_oligos::build(
            &prepared.plan,
            &prepared.target,
            &Default::default(),
        )
        .unwrap();
        assert!(oligos.complete);
        // Outer PCR can reuse branch primers; order records deduplicate DNA.
        assert_eq!(oligos.uses.len(), 6);
        assert!(
            !serde_json::to_string(&oligos.oligos)
                .unwrap()
                .contains("Seed stock")
        );
        native_request.origin = None;
        assert!(physical::preflight(&input, &native_request, Some(&frame)).is_err());
        let opened = Open {
            source: event.clone(),
            origin: Some(pin),
            epoch: 2,
            request,
            selection,
            check: checked.report().clone(),
            target: checked.target().clone(),
            prepared: Some(prepared),
        };
        let staged = stage(&frame, &opened).unwrap();
        let candidate = staged.candidate.as_ref().unwrap();
        assert_eq!(
            verification::verify(&candidate.plan, &opened.target).status,
            verification::Status::ConsistentUnderAssumptions
        );
        let saved = collagene_core::cloning_outputs::finalize(
            &staged.workspace,
            &candidate.plan,
            staged.target.clone(),
            staged.inputs.clone(),
            &collagene_core::cloning_outputs::Options {
                name: "Domesticated saved stock".into(),
                show_prepared: true,
                save_primers: true,
                allow_incomplete: false,
            },
        )
        .unwrap();
        assert!(saved.receipt.primer_attachment.is_some());
        let reopened =
            workspace::load_workspace(&workspace::save_workspace(&saved.workspace).unwrap())
                .unwrap();
        for original in &frame.documents {
            let retained = reopened
                .documents
                .iter()
                .find(|d| d.id == original.id)
                .unwrap();
            assert_eq!(retained.head_event, original.head_event);
            assert_eq!(retained.head, original.head);
        }
        let asset = reopened
            .attachments
            .values()
            .find(|a| a.kind == physical::REVIEW_EVIDENCE_KIND)
            .unwrap();
        let record: serde_json::Value = serde_json::from_slice(&asset.bytes).unwrap();
        assert_eq!(record["version"], 3);
        let restored_origin: MaterialOrigin =
            serde_json::from_value(record["origin"].clone()).unwrap();
        restored_origin
            .validate_for_workspace(&reopened, &event)
            .unwrap();
        assert_eq!(
            serde_json::to_vec(&restored_origin).unwrap(),
            serde_json::to_vec(&origin).unwrap()
        );
        let restored: physical::Prepared =
            serde_json::from_value(record["prepared"].clone()).unwrap();
        restored
            .validate_binding(&opened.prepared.as_ref().unwrap().context_revision)
            .unwrap();
        assert_eq!(
            restored.proposal_findings,
            opened.prepared.as_ref().unwrap().proposal_findings
        );
        let mut unpinned = opened.clone();
        unpinned.origin = None;
        assert!(stage(&frame, &unpinned).is_err());
        let mut alias_frame = frame.clone();
        let docs = DocumentHistory::restore_all(&frame).unwrap();
        let owner = docs.iter().find(|d| d.id == event.document_id).unwrap();
        let alias = owner
            .fork_event(&event.event_id, "equal-stock-other-event")
            .unwrap();
        let mut alias_open = opened.clone();
        alias_open.source = EventBinding {
            document_id: alias.id.clone(),
            event_id: alias.head_event().into(),
            revision: alias.head_revision().into(),
        };
        alias.contribute(&mut alias_frame).unwrap();
        assert!(stage(&alias_frame, &alias_open).is_err());
        let mut wrong = opened;
        wrong.target.sequence = collagene_core::Sequence::dna("AAAACCCC").unwrap();
        assert!(stage(&frame, &wrong).is_err());
    }
}
