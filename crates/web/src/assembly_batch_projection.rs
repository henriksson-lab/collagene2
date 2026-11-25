//! Borrowed projection of editable construct rows into the existing typed batch
//! contract. Invalid scientific text stays a row-local draft error.
use crate::{
    assembly_alternatives::{self as alternatives, Construct, Library},
    assembly_draft_model::{DigestGeometry, Draft, Slot},
    cloning_methods_model::{MethodsState, Mode, OverlapInput, OverlapPreparation},
    cloning_preparations,
};
use collagene_core::{
    assembly_junction_edit::Edit,
    combinatorial_batch::{self as batch, RowRequest},
    combinatorial_design::{
        self as design, ChoiceSelection, ConstructOverrides, Method, PoolSubstrate,
        PoolSubstrateKind, PreparationOverride, SelectedConstruct, TargetDeclaration,
    },
    run_lineage::EventBinding,
};
use serde::Serialize;
use std::io;

pub struct Projection {
    pub request: batch::Request,
    /// Raw declaration identity as well as typed request identity. Navigation is
    /// excluded, but changing invalid text or an uncaptured source invalidates it.
    pub revision: String,
}
fn source_label(name: &str, choice_id: &str) -> String {
    if !name.trim().is_empty() && name.len() <= 256 && !name.chars().any(char::is_control) {
        name.into()
    } else {
        // The exact source name remains in its immutable workspace revision.
        // A decorative typed-catalog label must not reject a valid long FASTA name.
        format!("Source alternative {choice_id}")
    }
}
fn selections(library: &Library, row: &Construct) -> Vec<ChoiceSelection> {
    library
        .bins
        .iter()
        .zip(&row.choices)
        .map(|(bin, id)| ChoiceSelection {
            bin_id: bin.id.clone(),
            choice_id: id.clone(),
        })
        .collect()
}
fn typed_row(
    library: &Library,
    row: &Construct,
    catalog: &design::Design,
) -> Result<SelectedConstruct, String> {
    for (bin, id) in library.bins.iter().zip(&row.choices) {
        if !catalog
            .bins
            .iter()
            .find(|entry| entry.id == bin.id)
            .is_some_and(|entry| entry.choices.iter().any(|choice| choice.id == *id))
        {
            return Err(format!(
                "Choose a captured source for bin {} before preparing this construct.",
                bin.label
            ));
        }
    }
    let state = &row.edits.methods;
    let mut preparation_state = state.clone();
    let method = match state.mode {
        Mode::Overlap => {
            let (overlaps, options) = state.profile()?;
            Method::Overlap { overlaps, options }
        }
        Mode::Ligation => Method::Ligation {
            allow_nicked: state.allow_nicked,
        },
        Mode::GoldenGate => {
            if let Some(pool) = &state.golden_gate_pool {
                if state.allow_nicked {
                    return Err("Finite Golden Gate pool selection requires sealed products; disable retained nicks for this construct.".into());
                }
                let options = pool.options(state.golden_gate.enzyme)?;
                let mut substrates = Vec::with_capacity(library.bins.len());
                for (index, (bin, raw)) in library
                    .bins
                    .iter()
                    .zip(&mut preparation_state.rows)
                    .enumerate()
                {
                    let choice = pool.choices.get(index).copied().unwrap_or(0);
                    let kind = match (pool.reaction, choice) {
                        (false, 0) | (true, 0) => PoolSubstrateKind::Prepared,
                        (true, 1) => {
                            if raw.treatment != crate::cloning_methods_model::Treatment::Keep {
                                return Err("Whole-original reaction substrates cannot silently discard a terminal treatment; choose Keep terminal chemistry.".into());
                            }
                            *raw = OverlapInput { document: raw.document, ..Default::default() };
                            PoolSubstrateKind::WholeSource
                        },
                        (true, 2) => {
                            if raw.treatment != crate::cloning_methods_model::Treatment::Keep {
                                return Err("Whole-amplicon reaction substrates require Keep terminal chemistry; post-digest terminal treatment does not apply to this substrate.".into());
                            }
                            if !matches!(raw.preparation, OverlapPreparation::ExactPcr | OverlapPreparation::PcrCatalogDigest) {
                                return Err("Choose actual PCR primers before using a whole amplicon.".into());
                            }
                            raw.preparation = OverlapPreparation::ExactPcr;
                            raw.fragment = "1".into();
                            raw.treatment = crate::cloning_methods_model::Treatment::Keep;
                            PoolSubstrateKind::Prepared
                        },
                        _ => return Err("Whole-original and whole-amplicon substrates require declared reaction mode.".into()),
                    };
                    substrates.push(PoolSubstrate {
                        bin_id: bin.id.clone(),
                        kind,
                    });
                }
                Method::GoldenGatePool {
                    options,
                    substrates,
                    selection: pool.candidate.as_ref().map(|intent| design::PoolSelection {
                        source_revision: intent.source_revision.clone(),
                        selection: collagene_core::prepared_golden_gate::Selection {
                            context_revision: intent.context_revision.clone(),
                            candidate_index: intent.candidate_index,
                            directed_sequence_key: intent.directed_sequence_key.clone(),
                            physical_sequence_key: intent.physical_sequence_key.clone(),
                        },
                    }),
                }
            } else {
                Method::GoldenGate {
                    enzyme: state.golden_gate.enzyme,
                    allow_nicked: state.allow_nicked,
                }
            }
        }
        Mode::Tailing => {
            return Err(
                "TA/GC tailing is not a source-bin assembly method; inspect this row separately."
                    .into(),
            );
        }
    };
    let (_, prepared) =
        cloning_preparations::prepared_rows(&preparation_state, library.bins.len())?;
    let preparations = library
        .bins
        .iter()
        .zip(prepared)
        .zip(&preparation_state.rows)
        .map(|((bin, fragment), raw)| PreparationOverride {
            bin_id: bin.id.clone(),
            preparation: Some(fragment.preparation),
            treatment: raw.treatment.operation(),
        })
        .collect();
    let value = SelectedConstruct {
        id: row.id.clone(),
        name: row.edits.name.clone(),
        choices: selections(library, row),
        target: TargetDeclaration::OriginalRegions {
            junction_edits: row.edits.junction_edits.clone(),
        },
        overrides: Some(ConstructOverrides {
            preparations,
            method: Some(method),
            circular: Some(row.edits.circular),
        }),
    };
    design::validate_construct(catalog, &value).map_err(|e| e.to_string())?;
    Ok(value)
}

#[derive(Serialize)]
struct SourceKey<'a> {
    label: &'a str,
    source: Option<&'a EventBinding>,
    start: usize,
    end: usize,
    reverse: bool,
    digest_geometry: Option<&'a DigestGeometry>,
}
fn source_key<'a>(slot: &'a Slot, label: &'a str) -> SourceKey<'a> {
    SourceKey {
        label,
        source: slot
            .captured
            .as_ref()
            .map(|capture| &capture.source.binding),
        start: slot.start,
        end: slot.end,
        reverse: slot.reverse,
        digest_geometry: slot.digest_geometry.as_ref(),
    }
}
#[derive(Serialize)]
struct ChoiceKey<'a> {
    id: &'a str,
    source: SourceKey<'a>,
    preparation: &'a OverlapInput,
}
#[derive(Serialize)]
struct BinKey<'a> {
    id: &'a str,
    label: &'a str,
    choices: Vec<ChoiceKey<'a>>,
}
#[derive(Serialize)]
struct RowKey<'a> {
    id: &'a str,
    choices: &'a [String],
    name: &'a str,
    circular: bool,
    methods: &'a MethodsState,
    junction_edits: &'a [Edit],
    sources: Vec<SourceKey<'a>>,
}
struct HashBytes {
    hash: blake3::Hasher,
    bytes: usize,
}
impl io::Write for HashBytes {
    fn write(&mut self, value: &[u8]) -> io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(value.len())
            .filter(|n| *n <= 32 * 1024 * 1024)
            .ok_or_else(|| io::Error::other("Assembly projection declaration exceeds32MiB"))?;
        self.hash.update(value);
        Ok(value.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
fn declaration_revision(stored: &Draft, request: &batch::Request) -> Result<String, String> {
    let library = stored
        .alternatives
        .as_ref()
        .ok_or("Enable source alternatives before preparing a construct batch.")?;
    let bins = library
        .bins
        .iter()
        .map(|bin| BinKey {
            id: &bin.id,
            label: &bin.label,
            choices: bin
                .choices
                .iter()
                .map(|choice| ChoiceKey {
                    id: &choice.id,
                    source: source_key(&choice.slot, &bin.label),
                    preparation: &choice.row,
                })
                .collect(),
        })
        .collect::<Vec<_>>();
    let rows = library
        .rows
        .iter()
        .map(|row| {
            let active = row.id == library.active_row;
            let sources = if active {
                stored
                    .slots
                    .iter()
                    .map(|slot| source_key(slot, &slot.label))
                    .collect()
            } else {
                library
                    .bins
                    .iter()
                    .zip(&row.choices)
                    .map(|(bin, id)| {
                        let choice = bin
                            .choices
                            .iter()
                            .find(|choice| choice.id == *id)
                            .ok_or("The construct source alternative is unavailable.")?;
                        Ok(source_key(&choice.slot, &bin.label))
                    })
                    .collect::<Result<Vec<_>, String>>()?
            };
            Ok(RowKey {
                id: &row.id,
                choices: &row.choices,
                name: if active {
                    &stored.name
                } else {
                    &row.edits.name
                },
                circular: if active {
                    stored.circular
                } else {
                    row.edits.circular
                },
                methods: if active {
                    &stored.methods
                } else {
                    &row.edits.methods
                },
                junction_edits: if active {
                    &stored.junction_edits
                } else {
                    &row.edits.junction_edits
                },
                sources,
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    let request_revision = batch::request_revision(request).map_err(|e| e.to_string())?;
    let mut output = HashBytes {
        hash: blake3::Hasher::new(),
        bytes: 0,
    };
    serde_json::to_writer(
        &mut output,
        &(
            "assembly-batch-projection-1",
            &library.id,
            bins,
            rows,
            request_revision,
        ),
    )
    .map_err(|e| e.to_string())?;
    Ok(output.hash.finalize().to_hex().to_string())
}

pub fn project(draft: &Draft) -> Result<Projection, String> {
    let stored = alternatives::for_storage(draft)?;
    let library = stored
        .alternatives
        .as_ref()
        .ok_or("Enable source alternatives before preparing a construct batch.")?;
    let pending = alternatives::require_active_choice_match(&stored).err();
    let bins = library
        .bins
        .iter()
        .map(|bin| design::Bin {
            id: bin.id.clone(),
            label: bin.label.clone(),
            choices: bin
                .choices
                .iter()
                .filter_map(|choice| {
                    let capture = choice.slot.captured.as_ref()?;
                    Some(design::Choice {
                        id: choice.id.clone(),
                        label: source_label(&capture.source.name, &choice.id),
                        source: capture.source.binding.clone(),
                        start: choice.slot.start,
                        end: choice.slot.end,
                        reverse: choice.slot.reverse,
                        preparation: None,
                        treatment: None,
                    })
                })
                .collect(),
        })
        .collect::<Vec<_>>();
    let catalog = if bins.iter().any(|bin| bin.choices.is_empty()) {
        None
    } else {
        Some(design::Design {
            id: library.id.clone(),
            name: format!("Assembly alternatives {}", library.id),
            circular: false,
            bins,
            method: Method::Ligation {
                allow_nicked: false,
            },
        })
    };
    if let Some(catalog) = &catalog {
        design::validate_design(catalog).map_err(|e| e.to_string())?;
    }
    let rows=library.rows.iter().map(|row|{
        let result=if row.id==library.active_row&&pending.is_some(){
            Err(pending.clone().expect("checked pending source"))
        }else if let Some(catalog)=&catalog{
            typed_row(library,row,catalog)
        }else{
            Err("Each ordered source bin needs a captured source before any construct can be prepared.".into())
        };
        match result{
            Ok(row)=>RowRequest::Ready(row),
            Err(message)=>RowRequest::DraftError{id:row.id.clone(),name:if row.id==library.active_row{stored.name.clone()}else{row.edits.name.clone()},choices:selections(library,row),message},
        }
    }).collect();
    let request = batch::Request {
        design: catalog,
        rows,
    };
    request.validate().map_err(|e| e.to_string())?;
    let revision = declaration_revision(&stored, &request)?;
    Ok(Projection { request, revision })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        assembly_draft_model::CapturedSource,
        cloning_methods_model::{OverlapPreparation, Treatment},
        golden_gate_pool_form::PoolForm,
        source_identity::SourceInfo,
    };
    use collagene_core::{
        Molecule, Topology,
        assembly_recipe::Preparation,
        chemistry::{EndTreatment, Terminus},
        workspace::Workspace,
    };
    use std::sync::Arc;

    fn source(name: &str, dna: &str) -> Arc<CapturedSource> {
        let molecule = Molecule::dna(name, dna, Topology::Linear).unwrap();
        let workspace = Workspace::from_documents(name, vec![molecule.clone()], 0).unwrap();
        let doc = &workspace.documents[0];
        Arc::new(CapturedSource {
            source: SourceInfo {
                binding: EventBinding {
                    document_id: doc.id.clone(),
                    event_id: doc.head_event.clone().unwrap(),
                    revision: doc.head.clone(),
                },
                name: name.into(),
                length: molecule.sequence.len(),
                alphabet: molecule.sequence.alphabet(),
                topology: molecule.topology,
            },
            molecule,
        })
    }
    fn ordinary() -> Draft {
        let mut draft = Draft::default();
        draft.name = "First construct".into();
        draft.circular = false;
        draft.methods.mode = Mode::Ligation;
        for (slot, captured) in draft
            .slots
            .iter_mut()
            .zip([source("backbone", "AACCGGTT"), source("insert", "TTGCAACC")])
        {
            slot.start = 0;
            slot.end = captured.molecule.sequence.len();
            slot.captured = Some(captured);
        }
        draft.sync_rows();
        for row in &mut draft.methods.rows {
            row.preparation = OverlapPreparation::DeclaredDuplex;
            row.left_terminus = Terminus::CLEAVED;
            row.right_terminus = Terminus::CLEAVED;
        }
        draft.validate().unwrap();
        draft
    }
    fn two_rows() -> Draft {
        let first = alternatives::enable(&ordinary(), "projection").unwrap();
        let selection = first.alternatives.as_ref().unwrap().rows[0].choices.clone();
        alternatives::add_row(&first, selection).unwrap()
    }
    fn ready(request: &batch::Request, index: usize) -> &SelectedConstruct {
        let RowRequest::Ready(row) = &request.rows[index] else {
            panic!("row is not ready")
        };
        row
    }

    #[test]
    fn every_construct_keeps_its_literal_preparations_method_topology_and_treatment() {
        let mut draft = two_rows();
        draft.name = "Second construct".into();
        draft.circular = true;
        draft.methods.mode = Mode::Overlap;
        draft.methods.overlaps = "4,4".into();
        let row = &mut draft.methods.rows[0];
        row.preparation = OverlapPreparation::ExactPcr;
        row.forward = "AACCGGTT".into();
        row.reverse_primer = "AACCGGTT".into();
        row.forward_tail = "TT".into();
        row.treatment = Treatment::Phosphorylate;
        draft.sync_rows();
        let result = project(&draft).unwrap();
        let first = ready(&result.request, 0);
        let second = ready(&result.request, 1);
        let first_settings = first.overrides.as_ref().unwrap();
        let second_settings = second.overrides.as_ref().unwrap();
        assert!(matches!(
            first_settings.method,
            Some(Method::Ligation {
                allow_nicked: false
            })
        ));
        assert_eq!(first_settings.circular, Some(false));
        assert!(matches!(
            &first_settings.preparations[0].preparation,
            Some(Preparation::DeclaredDuplex {
                left: Terminus::CLEAVED,
                right: Terminus::CLEAVED
            })
        ));
        assert!(
            matches!(&second_settings.method,Some(Method::Overlap{overlaps,..})if overlaps==&vec![4,4])
        );
        assert_eq!(second_settings.circular, Some(true));
        let Some(Preparation::ExactPcr { forward, reverse }) =
            &second_settings.preparations[0].preparation
        else {
            panic!("missing row-local PCR")
        };
        assert_eq!(forward.ordered_sequence(), "TTAACCGGTT");
        assert_eq!(reverse.ordered_sequence(), "AACCGGTT");
        assert_eq!(
            second_settings.preparations[0].treatment,
            Some(EndTreatment::Phosphorylate)
        );
        assert_ne!(first.id, second.id);
        assert_eq!(first.name, "First construct");
        assert_eq!(second.name, "Second construct");
    }

    #[test]
    fn unfinished_rows_do_not_hide_valid_rows_and_invalid_text_still_changes_the_key() {
        let mut draft = two_rows();
        draft.methods.rows[0].preparation = OverlapPreparation::ExactPcr;
        draft.methods.rows[0].forward = "?? unfinished".into();
        let broken = project(&draft).unwrap();
        assert!(matches!(&broken.request.rows[0], RowRequest::Ready(_)));
        assert!(matches!(
            &broken.request.rows[1],
            RowRequest::DraftError { .. }
        ));
        draft.methods.rows[0].preparation = OverlapPreparation::DeclaredDuplex;
        draft.methods.mode = Mode::Overlap;
        draft.methods.overlaps = "4".into();
        draft.methods.comparisons = "bad-one".into();
        let first = project(&draft).unwrap();
        draft.methods.comparisons = "bad-two".into();
        let second = project(&draft).unwrap();
        assert_eq!(
            batch::request_revision(&first.request).unwrap(),
            batch::request_revision(&second.request).unwrap()
        );
        assert_ne!(first.revision, second.revision);
        assert_eq!(draft.methods.comparisons, "bad-two");
    }

    #[test]
    fn navigation_keeps_declaration_identity_but_pending_active_sources_do_not_reuse_saved_choices()
    {
        let mut draft = two_rows();
        let original = project(&draft).unwrap();
        draft.active_slot = 1;
        assert_eq!(project(&draft).unwrap().revision, original.revision);
        let first_id = draft.alternatives.as_ref().unwrap().rows[0].id.clone();
        let other = alternatives::open_row(&draft, &first_id).unwrap();
        assert_eq!(project(&other).unwrap().revision, original.revision);
        let capture = source("pending replacement", "CCGTAAGG");
        draft.slots[1].captured = Some(capture);
        draft.slots[1].start = 0;
        draft.slots[1].end = 8;
        draft.sync_rows();
        let first = project(&draft).unwrap();
        assert!(matches!(&first.request.rows[0], RowRequest::Ready(_)));
        assert!(
            matches!(&first.request.rows[1],RowRequest::DraftError{message,..}if message.contains("capture"))
        );
        draft.slots[1].start = 1;
        draft.sync_rows();
        let second = project(&draft).unwrap();
        assert_eq!(
            batch::request_revision(&first.request).unwrap(),
            batch::request_revision(&second.request).unwrap()
        );
        assert_ne!(first.revision, second.revision);
        assert_eq!(draft.slots[1].start, 1);
    }

    #[test]
    fn unassigned_sources_fail_but_explicit_reaction_rows_remain_inspectable_without_a_pin() {
        let empty = alternatives::enable(&Draft::default(), "unassigned").unwrap();
        let result = project(&empty).unwrap();
        assert!(result.request.design.is_none());
        assert!(
            result
                .request
                .rows
                .iter()
                .all(|row| matches!(row, RowRequest::DraftError { .. }))
        );
        let mut draft = two_rows();
        draft.methods.mode = Mode::GoldenGate;
        draft.methods.golden_gate_pool = Some(PoolForm {
            reaction: true,
            choices: vec![1, 0],
            circular: false,
            ..PoolForm::default()
        });
        let result = project(&draft).unwrap();
        assert!(matches!(&result.request.rows[0], RowRequest::Ready(_)));
        let row = ready(&result.request, 1);
        let Method::GoldenGatePool {
            selection,
            substrates,
            options,
        } = row.overrides.as_ref().unwrap().method.as_ref().unwrap()
        else {
            panic!("expected declared reaction");
        };
        assert!(selection.is_none());
        assert!(!options.circular);
        assert!(matches!(substrates[0].kind, PoolSubstrateKind::WholeSource));
        assert!(matches!(substrates[1].kind, PoolSubstrateKind::Prepared));
        assert_eq!(
            draft.methods.golden_gate_pool.as_ref().unwrap().choices,
            vec![1, 0]
        );
        let long = "name".repeat(100);
        assert_eq!(
            source_label(&long, "choice-id"),
            "Source alternative choice-id"
        );
        assert_eq!(
            source_label("Original source", "choice-id"),
            "Original source"
        );
    }

    #[test]
    fn whole_amplicon_projection_keeps_primers_and_rejects_discarded_treatments_locally() {
        let mut draft = two_rows();
        draft.methods.mode = Mode::GoldenGate;
        draft.methods.golden_gate_pool = Some(PoolForm {
            reaction: true,
            choices: vec![2, 0],
            circular: false,
            ..PoolForm::default()
        });
        let raw = &mut draft.methods.rows[0];
        raw.preparation = OverlapPreparation::ExactPcr;
        raw.forward = "AACCGGTT".into();
        raw.reverse_primer = "AACCGGTT".into();
        let projection = project(&draft).unwrap();
        let row = ready(&projection.request, 1);
        let overrides = row.overrides.as_ref().unwrap();
        let Some(collagene_core::assembly_recipe::Preparation::ExactPcr { forward, reverse }) =
            &overrides.preparations[0].preparation
        else {
            panic!("whole amplicon did not retain PCR");
        };
        assert_eq!(forward.annealing.as_str(), "AACCGGTT");
        assert_eq!(reverse.annealing.as_str(), "AACCGGTT");
        draft.methods.rows[0].treatment = Treatment::Phosphorylate;
        let projection = project(&draft).unwrap();
        assert!(matches!(&projection.request.rows[0], RowRequest::Ready(_)));
        assert!(
            matches!(&projection.request.rows[1], RowRequest::DraftError { message, .. } if message.contains("terminal chemistry"))
        );
        draft.methods.rows[0].treatment = Treatment::Keep;
        draft.methods.allow_nicked = true;
        let projection = project(&draft).unwrap();
        assert!(
            matches!(&projection.request.rows[1], RowRequest::DraftError { message, .. } if message.contains("sealed products"))
        );
    }

    #[test]
    fn candidate_intent_roundtrips_per_construct_and_changes_declaration_identity() {
        use crate::golden_gate_pool_form::CandidateIntent;
        let mut draft = two_rows();
        draft.methods.mode = Mode::GoldenGate;
        draft.methods.golden_gate_pool = Some(PoolForm {
            reaction: true,
            choices: vec![1, 0],
            circular: false,
            ..PoolForm::default()
        });
        let before = project(&draft).unwrap();
        let intent = CandidateIntent {
            source_revision: "d".repeat(64),
            context_revision: "a".repeat(64),
            candidate_index: 3,
            directed_sequence_key: "b".repeat(64),
            physical_sequence_key: "c".repeat(64),
        };
        let saved = alternatives::set_active_pool_candidate(&draft, Some(intent.clone())).unwrap();
        let after = project(&saved).unwrap();
        assert_ne!(before.revision, after.revision);
        let Method::GoldenGatePool {
            selection: Some(selection),
            ..
        } = ready(&after.request, 1)
            .overrides
            .as_ref()
            .unwrap()
            .method
            .as_ref()
            .unwrap()
        else {
            panic!("missing explicit intent");
        };
        assert_eq!(selection.source_revision, intent.source_revision);
        assert_eq!(
            selection.selection.context_revision,
            intent.context_revision
        );
        assert_eq!(selection.selection.candidate_index, 3);
        let bytes = crate::assembly_draft_model::encode(&saved).unwrap();
        let restored = crate::assembly_draft_model::decode_version(&bytes, 3).unwrap();
        assert_eq!(project(&restored).unwrap().revision, after.revision);
        // Equal DNA under a different source event is a new declaration. The
        // old explicit pin remains old; projection must never refresh its token.
        let mut changed = saved.clone();
        let mut alias = changed.slots[0].captured.as_deref().unwrap().clone();
        alias.source.binding.event_id = format!("event-{}", "f".repeat(64));
        let alias = Arc::new(alias);
        changed.slots[0].captured = Some(alias.clone());
        let library = changed.alternatives.as_mut().unwrap();
        let active = library
            .rows
            .iter()
            .find(|row| row.id == library.active_row)
            .unwrap();
        let selected = active.choices[0].clone();
        library.bins[0]
            .choices
            .iter_mut()
            .find(|choice| choice.id == selected)
            .unwrap()
            .slot
            .captured = Some(alias);
        let aliased = project(&changed).unwrap();
        assert_ne!(aliased.revision, after.revision);
        let Method::GoldenGatePool {
            selection: Some(old_pin),
            ..
        } = ready(&aliased.request, 1)
            .overrides
            .as_ref()
            .unwrap()
            .method
            .as_ref()
            .unwrap()
        else {
            panic!("missing retained old pin");
        };
        assert_eq!(old_pin.source_revision, intent.source_revision);
        assert_eq!(old_pin.selection.context_revision, intent.context_revision);
        for old_version in [1, 2] {
            let mut value = serde_json::to_value(&saved).unwrap();
            value["version"] = serde_json::json!(old_version);
            assert!(
                crate::assembly_draft_model::decode(&serde_json::to_vec(&value).unwrap()).is_err()
            );
        }
        let cleared = alternatives::set_active_pool_candidate(&saved, None).unwrap();
        assert_eq!(project(&cleared).unwrap().revision, before.revision);
    }
}
