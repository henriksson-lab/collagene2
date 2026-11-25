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
        let expected = physical::revision(&input, &prepared.request).map_err(|e| e.to_string())?;
        prepared
            .validate_binding(&expected)
            .map_err(|e| e.to_string())?;
        if fingerprint(&prepared.request.design).map_err(|e| e.to_string())?
            != fingerprint(&opened.request).map_err(|e| e.to_string())?
            || fingerprint(&prepared.request.selection).map_err(|e| e.to_string())?
                != fingerprint(&opened.selection).map_err(|e| e.to_string())?
            || prepared.target != opened.target
            || prepared.preparation.plan.inputs.len() != 1
            || prepared.preparation.plan.inputs[0].id != physical::INPUT_ID
            || prepared.preparation.plan.inputs[0].molecule != source
        {
            return Err("The proposed physical preparation differs from the checked original source and independent intent.".into());
        }
        Some(Candidate::new(prepared.preparation.plan.clone()))
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
    let evidence=serde_json::to_vec(&serde_json::json!({"version":1,"source":opened.source,"check":checked.report(),"authority":"Checked reference constraints only; fresh reaction verification remains required."})).map_err(|e|e.to_string())?;
    if evidence.len() > collagene_core::workspace::MAX_ENTRY_BYTES {
        return Err("Domestication evidence exceeds its attachment bound.".into());
    }
    let hash = blake3::hash(&evidence).to_hex().to_string();
    staged
        .attachments
        .entry(hash)
        .or_insert(WorkspaceAttachment {
            kind: "domestication-intent-v1".into(),
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
}
