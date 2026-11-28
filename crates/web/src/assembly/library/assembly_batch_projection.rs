//! Borrowed projection of editable construct rows into the existing typed batch
//! contract. Invalid scientific text stays a row-local draft error.
use crate::assembly::draft::assembly_alternatives as alternatives;
use crate::assembly::draft::assembly_alternatives::Construct;
use crate::assembly::draft::assembly_alternatives::Library;
use crate::assembly::draft::model::DigestGeometry;
use crate::assembly::draft::model::Draft;
use crate::assembly::draft::model::Slot;
use crate::cloning::preparation::cloning_methods_model::MethodsState;
use crate::cloning::preparation::cloning_methods_model::Mode;
use crate::cloning::preparation::cloning_methods_model::OverlapInput;
use crate::cloning::preparation::cloning_methods_model::OverlapPreparation;
use crate::cloning::preparation::cloning_preparations;
use collagene_core::{
    cloning::assembly::junction_edit::Edit,
    material::lineage::EventBinding,
    planning::combinatorial::{
        self as design, ChoiceSelection, ConstructOverrides, Method, PoolSubstrate,
        PoolSubstrateKind, PreparationOverride, SelectedConstruct, TargetDeclaration,
        batch::{self, RowRequest},
    },
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
            let finishing = state.overlap_finishing(overlaps.clone(), state.rows.len())?;
            Method::Overlap {
                overlaps,
                options,
                finishing,
            }
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
                            if raw.treatment != crate::cloning::preparation::cloning_methods_model::Treatment::Keep {
                                return Err("Whole-original reaction substrates cannot silently discard a terminal treatment; choose Keep terminal chemistry.".into());
                            }
                            *raw = OverlapInput { document: raw.document, run_origin:raw.run_origin.clone(), source_generation:raw.source_generation.clone(), ..Default::default() };
                            PoolSubstrateKind::WholeSource
                        },
                        (true, 2) => {
                            if raw.treatment != crate::cloning::preparation::cloning_methods_model::Treatment::Keep {
                                return Err("Whole-amplicon reaction substrates require Keep terminal chemistry; post-digest terminal treatment does not apply to this substrate.".into());
                            }
                            raw.use_complete_pcr_amplicon()?;
                            raw.treatment = crate::cloning::preparation::cloning_methods_model::Treatment::Keep;
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
                        selection: collagene_core::cloning::golden_gate::prepared::Selection {
                            context_revision: intent.context_revision.clone(),
                            candidate_index: intent.candidate_index,
                            directed_sequence_key: intent.directed_sequence_key.clone(),
                            physical_sequence_key: intent.physical_sequence_key.clone(),
                        },
                    }),
                }
            } else {
                Method::OrderedTypeIis {
                    enzyme: state.golden_gate.enzyme,
                    allow_nicked: state.allow_nicked,
                }
            }
        }
        Mode::MultiSiteReference => {
            if !(2..=5).contains(&library.bins.len()) {
                return Err("Reference exchange requires2..5 whole original parents.".into());
            }
            for raw in &state.rows {
                if raw.preparation != OverlapPreparation::WholeLinear
                    || raw.treatment
                        != crate::cloning::preparation::cloning_methods_model::Treatment::Keep
                    || raw.repair.is_some()
                {
                    return Err("Explicitly choose whole original parents without preparation, repair or treatment for reference exchange.".into());
                }
            }
            let form = state
                .multisite
                .as_ref()
                .ok_or("Configure the reference profile first.")?;
            Method::MultiSiteReference {
                request: form.request()?,
                selection: form.selection(),
            }
        }
        Mode::Tailing => {
            if library.bins.len() != 2 || state.rows.len() != 2 {
                return Err(
                    "TA/GC library constructs require two explicit prepared source bins.".into(),
                );
            }
            Method::Tailing {
                profile: state.tailing_profile,
                actions: library
                    .bins
                    .iter()
                    .enumerate()
                    .map(|(index, bin)| design::TailingAction {
                        bin_id: bin.id.clone(),
                        action: if index == 0 {
                            state.vector.action
                        } else {
                            state.insert.action
                        },
                    })
                    .collect(),
                allow_nicked: state.allow_nicked,
            }
        }
    };
    let ordered_preparations = !matches!(
        &method,
        Method::GoldenGatePool { .. } | Method::MultiSiteReference { .. }
    );
    let preparations = if state.mode == Mode::MultiSiteReference {
        library
            .bins
            .iter()
            .zip(&state.rows)
            .map(|(bin, raw)| {
                Ok(PreparationOverride {
                    bin_id: bin.id.clone(),
                    preparation_reverse: ordered_preparations && raw.reverse,
                    preparation: None,
                    treatment: None,
                    generation: raw.generation()?,
                })
            })
            .collect::<Result<Vec<_>, String>>()?
    } else {
        let (_, prepared) = cloning_preparations::prepared_rows(
            &preparation_state,
            preparation_state.all_input_ids()?.len(),
        )?;
        library
            .bins
            .iter()
            .zip(prepared)
            .zip(&preparation_state.rows)
            .map(|((bin, fragment), raw)| PreparationOverride {
                bin_id: bin.id.clone(),
                preparation_reverse: ordered_preparations && fragment.reverse,
                generation: fragment.generation,
                preparation: Some(fragment.preparation),
                treatment: raw.treatment.operation(),
            })
            .collect()
    };
    let value = SelectedConstruct {
        id: row.id.clone(),
        name: row.edits.name.clone(),
        choices: selections(library, row),
        target: TargetDeclaration::OriginalRegions {
            junction_edits: row.edits.junction_edits.clone(),
            substitutions: row.edits.substitutions.clone(),
            edits: row.edits.edits.clone(),
        },
        overrides: Some(ConstructOverrides {
            preparations,
            method: Some(method),
            circular: Some(row.edits.circular),
        }),
    };
    catalog
        .validate_construct(&value)
        .map_err(|e| e.to_string())?;
    Ok(value)
}

#[derive(Serialize)]
struct SourceKey<'a> {
    label: &'a str,
    source: Option<&'a EventBinding>,
    anneal_bottom: Option<&'a EventBinding>,
    selector: &'a collagene_core::cloning::assembly::intent::Selector,
    material_origin: Option<&'a str>,
    digest_geometry: Option<&'a DigestGeometry>,
}
fn source_key<'a>(slot: &'a Slot, label: &'a str) -> SourceKey<'a> {
    SourceKey {
        label,
        source: slot
            .captured
            .as_ref()
            .map(|capture| &capture.source.binding),
        anneal_bottom: slot
            .anneal_bottom
            .as_ref()
            .map(|capture| &capture.source.binding),
        selector: &slot.selector,
        material_origin: slot.material_origin.as_ref().map(|pin| pin.revision()),
        digest_geometry: slot.digest_geometry.as_ref(),
    }
}
#[derive(Serialize)]
struct ChoiceKey<'a> {
    id: &'a str,
    source: SourceKey<'a>,
    preparation: OverlapInput,
    physical_origin: Option<&'a str>,
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
    methods: MethodsState,
    physical_origins: Vec<Option<&'a str>>,
    junction_edits: &'a [Edit],
    substitutions: &'a [collagene_core::cloning::assembly::substitution_intent::Substitution],
    edits: &'a [collagene_core::cloning::assembly::edit_intent::Edit],
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
/// Raw declaration identity for UI guards. Large immutable captures are represented
/// by validated source revisions and cached origin identities, not copied/replayed.
pub fn declaration_revision(stored: &Draft) -> Result<String, String> {
    let library = stored
        .alternatives
        .as_ref()
        .ok_or("Enable source alternatives before preparing a construct batch.")?;
    if stored.slots.len() > 12
        || library.bins.is_empty()
        || library.bins.len() > alternatives::MAX_BINS
        || library.rows.len() > alternatives::MAX_ROWS
        || library
            .bins
            .iter()
            .any(|bin| bin.choices.len() > alternatives::MAX_CHOICES_PER_BIN)
        || library
            .bins
            .iter()
            .map(|bin| bin.choices.len())
            .sum::<usize>()
            > alternatives::MAX_CHOICES
    {
        return Err("Library declaration exceeds bounded row/source counts".into());
    }
    if stored.methods.rows.len() > 12
        || library
            .rows
            .iter()
            .any(|row| row.edits.methods.rows.len() > 12)
    {
        return Err("Library method rows exceed the declaration bound".into());
    }
    crate::assembly::draft::model::validate_controls(
        &stored.name,
        stored.circular,
        stored.active_slot,
        stored.slots.len(),
        &stored.methods,
        &stored.junction_edits,
    )?;
    for row in &library.rows {
        crate::assembly::draft::model::validate_controls(
            &row.edits.name,
            row.edits.circular,
            row.edits.active_slot,
            library.bins.len(),
            &row.edits.methods,
            &row.edits.junction_edits,
        )?;
    }
    for choice in library.bins.iter().flat_map(|bin| &bin.choices) {
        crate::assembly::draft::model::validate_preparation_input(&choice.row)?;
    }
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
                    preparation: {
                        let mut fields = choice.row.clone();
                        fields.run_origin = None;
                        fields
                    },
                    physical_origin: choice.row.run_origin.as_ref().map(|pin| pin.revision()),
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
                methods: {
                    let mut fields = if active {
                        stored.methods.clone()
                    } else {
                        row.edits.methods.clone()
                    };
                    for raw in &mut fields.rows {
                        raw.run_origin = None;
                    }
                    fields
                },
                physical_origins: if active {
                    &stored.methods
                } else {
                    &row.edits.methods
                }
                .rows
                .iter()
                .map(|raw| raw.run_origin.as_ref().map(|pin| pin.revision()))
                .collect(),
                junction_edits: if active {
                    &stored.junction_edits
                } else {
                    &row.edits.junction_edits
                },
                substitutions: if active {
                    &stored.substitutions
                } else {
                    &row.edits.substitutions
                },
                edits: if active {
                    &stored.edits
                } else {
                    &row.edits.edits
                },
                sources,
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    let mut output = HashBytes {
        hash: blake3::Hasher::new(),
        bytes: 0,
    };
    serde_json::to_writer(
        &mut output,
        &(
            "assembly-batch-projection-3",
            crate::assembly::draft::model::CURRENT_VERSION,
            collagene_core::verification::ENGINE_VERSION,
            &library.id,
            bins,
            rows,
        ),
    )
    .map_err(|e| e.to_string())?;
    Ok(output.hash.finalize().to_hex().to_string())
}

fn material_bindings(library: &Library) -> Result<Vec<design::LibraryMaterialBinding>, String> {
    let mut sources: Vec<(
        &EventBinding,
        &crate::history::run_origin_model::PinnedOrigin,
    )> = Vec::new();
    for choice in library.bins.iter().flat_map(|bin| &bin.choices) {
        let Some(origin) = choice.slot.material_origin.as_ref() else {
            continue;
        };
        let captured = choice
            .slot
            .captured
            .as_ref()
            .ok_or("Saved material source has no captured event")?;
        if origin.event() != &captured.source.binding {
            return Err("Library source material pin belongs to another exact event".into());
        }
        if let Some((_, previous)) = sources
            .iter()
            .find(|(event, _)| *event == &captured.source.binding)
        {
            if *previous != origin {
                return Err("Library source event has conflicting material origins".into());
            }
        } else {
            sources.push((&captured.source.binding, origin));
        }
    }
    collagene_core::material::capsule::material_capsules_preflight(sources.iter().map(
        |(_, origin)| match origin.as_ref() {
            collagene_core::material::run::MaterialOrigin::Run(c) => {
                collagene_core::material::capsule::MaterialOriginRef::Run(c)
            }
            collagene_core::material::run::MaterialOrigin::Derived(c) => {
                collagene_core::material::capsule::MaterialOriginRef::Derived(c)
            }
        },
    ))
    .map_err(|e| e.to_string())?;
    Ok(sources
        .into_iter()
        .map(|(event, origin)| design::LibraryMaterialBinding {
            event: event.clone(),
            origin: origin.as_ref().clone(),
        })
        .collect())
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
                        anneal_bottom: choice
                            .slot
                            .anneal_bottom
                            .as_ref()
                            .map(|capture| capture.source.binding.clone()),
                        selector: choice.slot.selector.clone(),
                        preparation_reverse: choice.row.reverse,
                        preparation: None,
                        treatment: None,
                        generation: None,
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
            material_bindings: material_bindings(library)?,
            name: format!("Assembly alternatives {}", library.id),
            circular: false,
            bins,
            method: Method::Ligation {
                allow_nicked: false,
            },
        })
    };
    if let Some(catalog) = &catalog {
        catalog.validate().map_err(|e| e.to_string())?;
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
    let revision = declaration_revision(draft)?;
    Ok(Projection { request, revision })
}

/// Same typed declaration for the ordinary dock and an opened library construct.
/// The ephemeral one-choice catalog has namespaced stable identities, never DNA-derived aliases.
pub fn active(draft: &Draft) -> Result<(design::Design, SelectedConstruct), String> {
    let ordinary;
    let current = if draft.alternatives.is_none() {
        ordinary =
            alternatives::enable(draft, crate::cloning::methods::multisite_form::ORDINARY_ID)?;
        &ordinary
    } else {
        draft
    };
    let projection = project(current)?;
    let catalog = projection
        .request
        .design
        .ok_or("Capture all original sources first.")?;
    let id = &current
        .alternatives
        .as_ref()
        .ok_or("Active construct unavailable.")?
        .active_row;
    let row = projection
        .request
        .rows
        .into_iter()
        .find_map(|row| match row {
            RowRequest::Ready(row) if &row.id == id => Some(Ok(row)),
            RowRequest::DraftError {
                id: row_id,
                message,
                ..
            } if &row_id == id => Some(Err(message)),
            _ => None,
        })
        .ok_or("Active construct is unavailable.")??;
    Ok((catalog, row))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::draft::model::CapturedSource;
    use crate::cloning::preparation::cloning_methods_model::OverlapPreparation;
    use crate::cloning::preparation::cloning_methods_model::Treatment;
    use crate::cloning::preparation::golden_gate_pool_form::PoolForm;
    use crate::sources::identity::SourceInfo;
    use collagene_core::{
        chemistry::{EndTreatment, Terminus},
        cloning::assembly::recipe::Preparation,
        documents::workspace::Workspace,
    };
    use collagene_model::{Molecule, Topology};
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
            {
                let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                    start: value,
                    ..
                } = &mut slot.selector
                else {
                    panic!("fixture requires reference selector")
                };
                *value = 0;
            }
            {
                let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                    end: value,
                    ..
                } = &mut slot.selector
                else {
                    panic!("fixture requires reference selector")
                };
                *value = captured.molecule.sequence.len();
            }
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
    fn requested_processing_is_a_row_local_error_and_raw_policy_is_retained() {
        for three_prime in [false, true] {
            let mut draft = two_rows();
            draft.methods.mode = Mode::Overlap;
            if three_prime {
                let mut form =
                    crate::cloning::preparation::cloning_methods_model::AnnealingForm::default();
                form.enabled = true;
                form.cuts = "unfinished 3-prime cuts".into();
                draft.methods.annealing = Some(form);
            } else {
                let mut form =
                    crate::cloning::preparation::cloning_methods_model::HomologyForm::default();
                form.enabled = true;
                form.exact_homology = false;
                form.cuts = "unfinished 5-prime cuts".into();
                draft.methods.homology = Some(form);
            }
            draft.sync_rows();
            let before = serde_json::to_vec(&draft).unwrap();
            let projection = project(&draft).unwrap();
            assert!(matches!(&projection.request.rows[0], RowRequest::Ready(_)));
            assert!(
                matches!(&projection.request.rows[1],RowRequest::DraftError{message,..} if !message.is_empty())
            );
            assert_eq!(serde_json::to_vec(&draft).unwrap(), before);
            if three_prime {
                draft.methods.annealing.as_mut().unwrap().enabled = false;
            } else {
                draft.methods.homology.as_mut().unwrap().enabled = false;
            }
            assert!(
                matches!(
                    &project(&draft).unwrap().request.rows[1],
                    RowRequest::Ready(_)
                ),
                "Inactive raw settings must remain saveable without requesting processing"
            );
        }
    }
    #[test]
    fn explicit_annealing_policy_and_planned_rescue_survive_row_projection_and_draft_recovery() {
        let mut draft = two_rows();
        draft.methods.mode = Mode::Overlap;
        draft.methods.method = collagene_core::cloning::overlap::profiles::OverlapMethod::InFusion;
        draft.methods.homology = None;
        let mut form = crate::cloning::preparation::cloning_methods_model::AnnealingForm::default();
        form.enabled = true;
        form.cuts = "15:15, 15:15".into();
        form.min_paired_bases = "4".into();
        form.planned_ecoli_rescue = true;
        draft.methods.annealing = Some(form);
        draft.sync_rows();
        let before = project(&draft).unwrap();
        let Some(Method::Overlap {
            finishing:
                collagene_core::cloning::assembly::recipe::OverlapFinishing::HomologyAnnealing {
                    policy,
                },
            ..
        }) = &ready(&before.request, 1).overrides.as_ref().unwrap().method
        else {
            panic!("actual annealing policy must be retained");
        };
        assert_eq!(policy.cuts[0].left, 15);
        assert_eq!(
            policy.planned_rescue,
            Some(collagene_core::cloning::overlap::annealing::PlannedRescue::EcoliTransformation)
        );
        let restored = crate::assembly::draft::model::decode(
            &crate::assembly::draft::model::encode(&draft).unwrap(),
        )
        .unwrap();
        assert_eq!(project(&restored).unwrap().revision, before.revision);
        draft
            .methods
            .annealing
            .as_mut()
            .unwrap()
            .planned_ecoli_rescue = false;
        let invalid = project(&draft).unwrap();
        assert_ne!(invalid.revision, before.revision);
        assert!(matches!(
            &invalid.request.rows[1],
            RowRequest::DraftError { .. }
        ));
        assert!(matches!(&invalid.request.rows[0], RowRequest::Ready(_)));
    }

    #[test]
    fn row_local_literals_enter_the_native_declaration_and_raw_projection_key() {
        let draft = two_rows();
        let draft =
            crate::assembly::edit::assembly_edit_model::declare(&draft, 0, "3", "5", "TAG", None)
                .unwrap();
        let draft =
            crate::assembly::edit::assembly_substitution_model::declare(&draft, 1, "1", "A", None)
                .unwrap();
        let projected = project(&draft).unwrap();
        assert_eq!(declaration_revision(&draft).unwrap(), projected.revision);
        let TargetDeclaration::OriginalRegions {
            substitutions: first_subs,
            edits: first_edits,
            ..
        } = &ready(&projected.request, 0).target
        else {
            panic!("literal declaration")
        };
        assert!(first_subs.is_empty() && first_edits.is_empty());
        let TargetDeclaration::OriginalRegions {
            substitutions,
            edits,
            ..
        } = &ready(&projected.request, 1).target
        else {
            panic!("literal declaration")
        };
        assert_eq!(edits[0].expected, "CC");
        assert_eq!(edits[0].replacement, "TAG");
        assert_eq!(substitutions[0].expected, "T");
        assert_eq!(substitutions[0].replacement, "A");
        let mut changed = draft.clone();
        changed.edits[0].expected = "AA".into();
        let stale = project(&changed).unwrap();
        assert_ne!(projected.revision, stale.revision);
        assert_eq!(
            serde_json::to_vec(&ready(&projected.request, 0)).unwrap(),
            serde_json::to_vec(&ready(&stale.request, 0)).unwrap()
        );
        let mut changed = draft.clone();
        changed.edits[0].replacement = "unfinished?".into();
        assert_ne!(declaration_revision(&changed).unwrap(), projected.revision);
        assert!(project(&changed).is_err());
    }
    #[test]
    fn resolved_generation_is_explicit_in_overrides_and_partial_text_stays_row_local() {
        use crate::cloning::preparation::cloning_methods_model::SourceGenerationInput;
        use collagene_core::chemistry::methylation::MethylationState;
        let mut original = ordinary();
        for row in &mut original.methods.rows {
            row.preparation = OverlapPreparation::WholeLinear;
        }
        original.methods.rows[0].source_generation = Some(SourceGenerationInput {
            host: None,
            dam: "fully".into(),
            dcm: "unmodified".into(),
            cpg: "unknown".into(),
        });
        let library = alternatives::enable(&original, "genproj").unwrap();
        let projection = project(&library).unwrap();
        assert!(
            projection
                .request
                .design
                .as_ref()
                .unwrap()
                .bins
                .iter()
                .flat_map(|b| &b.choices)
                .all(|c| c.generation.is_none())
        );
        let overrides = &ready(&projection.request, 0)
            .overrides
            .as_ref()
            .unwrap()
            .preparations;
        assert!(overrides[0].generation.as_ref().unwrap().dam == MethylationState::FullyModified);
        assert!(overrides[1].generation.is_none());
        let mut cleared = library.clone();
        cleared.methods.rows[0].source_generation = None;
        let projection = project(&cleared).unwrap();
        assert!(
            ready(&projection.request, 0)
                .overrides
                .as_ref()
                .unwrap()
                .preparations[0]
                .generation
                .is_none()
        );
        let mut unfinished = library;
        unfinished.methods.rows[0]
            .source_generation
            .as_mut()
            .unwrap()
            .dam = "partial".into();
        unfinished.validate().unwrap();
        let projection = project(&unfinished).unwrap();
        assert!(
            matches!(&projection.request.rows[0],RowRequest::DraftError{message,..} if message.contains("source-generation"))
        );
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
        let Some(Preparation::ExactPcr {
            template_strand: collagene_core::chemistry::Strand::Top,
            forward,
            reverse,
        }) = &second_settings.preparations[0].preparation
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
            *value = 8;
        }
        draft.sync_rows();
        let first = project(&draft).unwrap();
        assert!(matches!(&first.request.rows[0], RowRequest::Ready(_)));
        assert!(
            matches!(&first.request.rows[1],RowRequest::DraftError{message,..}if message.contains("capture"))
        );
        {
            let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                start: value,
                ..
            } = &mut draft.slots[1].selector
            else {
                panic!("fixture requires reference selector")
            };
            *value = 1;
        }
        draft.sync_rows();
        let second = project(&draft).unwrap();
        assert_eq!(
            batch::request_revision(&first.request).unwrap(),
            batch::request_revision(&second.request).unwrap()
        );
        assert_ne!(first.revision, second.revision);
        assert_eq!(draft.slots[1].selector.reference_interval().unwrap().0, 1);
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
        let Some(collagene_core::cloning::assembly::recipe::Preparation::ExactPcr {
            template_strand: collagene_core::chemistry::Strand::Top,
            forward,
            reverse,
        }) = &overrides.preparations[0].preparation
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
    fn whole_amplicon_projection_preserves_intentional_mismatch_policy() {
        let mut draft = two_rows();
        draft.methods.mode = Mode::GoldenGate;
        draft.methods.golden_gate_pool = Some(PoolForm {
            reaction: true,
            choices: vec![2, 0],
            circular: false,
            ..PoolForm::default()
        });
        let raw = &mut draft.methods.rows[0];
        raw.preparation = OverlapPreparation::MismatchPcr;
        raw.forward = "ATGA".into();
        raw.reverse_primer = "GTCA".into();
        raw.mismatch_max_mismatches = "1".into();
        raw.mismatch_three_prime_exact = "2".into();
        let before = draft.clone();
        let projection = project(&draft).unwrap();
        let row = ready(&projection.request, 1);
        let overrides = row.overrides.as_ref().unwrap();
        let Some(collagene_core::cloning::assembly::recipe::Preparation::MismatchPcr {
            forward,
            reverse,
            options,
        }) = &overrides.preparations[0].preparation
        else {
            panic!("whole amplicon discarded intentional mismatches");
        };
        assert_eq!(forward.annealing.as_str(), "ATGA");
        assert_eq!(reverse.annealing.as_str(), "GTCA");
        assert_eq!(options.binding.max_mismatches, 1);
        assert_eq!(options.binding.min_three_prime_exact, 2);
        assert!(
            draft == before,
            "projection must not rewrite the declared target or draft"
        );
        draft.methods.rows[0].mismatch_max_mismatches = "2e".into();
        let unfinished = project(&draft).unwrap();
        assert_ne!(unfinished.revision, projection.revision);
        assert!(matches!(unfinished.request.rows[0], RowRequest::Ready(_)));
        assert!(matches!(
            unfinished.request.rows[1],
            RowRequest::DraftError { .. }
        ));
    }

    #[test]
    fn candidate_intent_roundtrips_per_construct_and_changes_declaration_identity() {
        use crate::cloning::preparation::golden_gate_pool_form::CandidateIntent;
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
        let bytes = crate::assembly::draft::model::encode(&saved).unwrap();
        let restored = crate::assembly::draft::model::decode(&bytes).unwrap();
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
                crate::assembly::draft::model::decode(&serde_json::to_vec(&value).unwrap())
                    .is_err()
            );
        }
        let cleared = alternatives::set_active_pool_candidate(&saved, None).unwrap();
        assert_eq!(project(&cleared).unwrap().revision, before.revision);
    }
    fn reference_draft() -> Draft {
        let mut draft = ordinary();
        draft.methods.mode = Mode::MultiSiteReference;
        let mut form = crate::cloning::methods::multisite_form::Form::default();
        form.reviewed = true;
        draft.methods.multisite = Some(form);
        for row in &mut draft.methods.rows {
            row.preparation = OverlapPreparation::WholeLinear;
            row.treatment = Treatment::Keep;
            row.repair = None;
        }
        draft.sync_rows();
        draft
    }
    #[test]
    fn ordinary_reference_projection_has_matching_roles_and_explicit_no_preparation_without_a_saved_matrix()
     {
        let draft = reference_draft();
        let before = serde_json::to_vec(&draft).unwrap();
        let (catalog, row) = active(&draft).unwrap();
        assert!(draft.alternatives.is_none());
        assert_eq!(before, serde_json::to_vec(&draft).unwrap());
        let roles = draft.methods.role_ids(draft.slots.len()).unwrap();
        for (i, bin) in catalog.bins.iter().enumerate() {
            let wanted =
                design::source_role_ids(&catalog.id, &row.id, &bin.id, &row.choices[i].choice_id)
                    .unwrap();
            assert_eq!((&roles[i].input_id, &roles[i].role), (&wanted.0, &wanted.1));
            assert_eq!(
                bin.choices[0].source,
                draft.slots[i].captured.as_ref().unwrap().source.binding
            );
            assert_eq!(
                (
                    bin.choices[0].selector.reference_interval().unwrap().0,
                    bin.choices[0].selector.reference_interval().unwrap().1,
                    bin.choices[0].selector.reference_interval().unwrap().2
                ),
                (
                    draft.slots[i].selector.reference_interval().unwrap().0,
                    draft.slots[i].selector.reference_interval().unwrap().1,
                    draft.slots[i].selector.reference_interval().unwrap().2
                )
            );
        }
        let overrides = row.overrides.unwrap();
        assert!(matches!(
            overrides.method,
            Some(Method::MultiSiteReference {
                selection: None,
                ..
            })
        ));
        assert!(
            overrides
                .preparations
                .iter()
                .all(|p| p.preparation.is_none() && p.treatment.is_none())
        );
    }
    #[test]
    fn reference_projection_refuses_to_discard_an_existing_preparation_or_treatment() {
        let mut draft = reference_draft();
        draft.methods.rows[0].preparation = OverlapPreparation::ExactPcr;
        draft.methods.rows[0].forward = "ACGT".into();
        assert!(active(&draft).unwrap_err().contains("without preparation"));
        assert_eq!(draft.methods.rows[0].forward, "ACGT");
        draft.methods.rows[0].preparation = OverlapPreparation::WholeLinear;
        draft.methods.rows[0].treatment = Treatment::Phosphorylate;
        assert!(active(&draft).unwrap_err().contains("without preparation"));
        draft.methods.rows[0].treatment = Treatment::Keep;
        assert!(active(&draft).is_ok());
    }
    #[test]
    fn reference_library_pin_and_unfinished_profile_persist_as_intent_without_proof() {
        let mut draft = alternatives::enable(&reference_draft(), "ref-library").unwrap();
        let form = draft.methods.multisite.as_mut().unwrap();
        form.selection = Some(crate::cloning::methods::multisite_form::ComponentIntent {
            batch_revision: "a".repeat(64),
            candidate_id: 2,
            component_id: 1,
            source_revision: "b".repeat(64),
        });
        let saved = alternatives::save_active(&draft).unwrap();
        let bytes = crate::assembly::draft::model::encode_raw(&saved).unwrap();
        let restored = crate::assembly::draft::model::decode(&bytes).unwrap();
        let (_, row) = active(&restored).unwrap();
        let Some(Method::MultiSiteReference {
            selection: Some(pin),
            ..
        }) = row.overrides.unwrap().method
        else {
            panic!("reference component intent lost")
        };
        assert_eq!(pin.source_revision, "b".repeat(64));
        assert_eq!(pin.selection.candidate_id, 2);
        let mut raw = restored;
        raw.methods.multisite.as_mut().unwrap().profile_json = "{unfinished profile Ω".into();
        let raw = alternatives::save_active(&raw).unwrap();
        let restored = crate::assembly::draft::model::decode(
            &crate::assembly::draft::model::encode_raw(&raw).unwrap(),
        )
        .unwrap();
        assert_eq!(
            restored.methods.multisite.as_ref().unwrap().profile_json,
            "{unfinished profile Ω"
        );
        assert!(
            active(&restored).is_err(),
            "editable malformed profile is retained but cannot generate a scientific request"
        );
    }
    #[test]
    fn ordinary_reference_and_edited_target_metadata_equal_the_frozen_intent_builder() {
        use collagene_core::{
            cloning::assembly::{
                intent as assembly_intent, junction_edit as assembly_junction_edit,
            },
            planning::plan::PlanInput,
        };
        let originals = vec![
            Molecule::dna("first source", "AACCGGTTAA", Topology::Linear).unwrap(),
            Molecule::dna("second source", "GGTTAACCGG", Topology::Linear).unwrap(),
        ];
        let frame = Workspace::from_documents("Real source owners", originals.clone(), 0).unwrap();
        let mut draft = reference_draft();
        for (i, slot) in draft.slots.iter_mut().enumerate() {
            let doc = &frame.documents[i];
            slot.captured = Some(Arc::new(CapturedSource {
                source: SourceInfo {
                    binding: EventBinding {
                        document_id: doc.id.clone(),
                        event_id: doc.head_event.clone().unwrap(),
                        revision: doc.head.clone(),
                    },
                    name: originals[i].name.clone(),
                    length: originals[i].sequence.len(),
                    alphabet: originals[i].sequence.alphabet(),
                    topology: originals[i].topology,
                },
                molecule: originals[i].clone(),
            }));
            {
                let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                    start: value,
                    ..
                } = &mut slot.selector
                else {
                    panic!("fixture requires reference selector")
                };
                *value = i + 1;
            }
            {
                let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                    end: value,
                    ..
                } = &mut slot.selector
                else {
                    panic!("fixture requires reference selector")
                };
                *value = 8 + i;
            }
            {
                let collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                    reverse: value,
                    ..
                } = &mut slot.selector
                else {
                    panic!("fixture requires reference selector")
                };
                *value = i == 1;
            }
        }
        draft.sync_rows();
        for edited in [false, true] {
            draft.junction_edits = if edited {
                vec![Edit {
                    junction: 0,
                    remove_left: 1,
                    remove_right: 1,
                    inserted: "GAT".into(),
                }]
            } else {
                vec![]
            };
            let roles = draft.methods.role_ids(draft.slots.len()).unwrap();
            let inputs = roles
                .iter()
                .zip(&draft.slots)
                .map(|(r, s)| {
                    let c = s.captured.as_ref().unwrap();
                    PlanInput {
                        id: r.input_id.clone(),
                        revision: c.source.binding.revision.clone(),
                        molecule: c.molecule.clone(),
                    }
                })
                .collect::<Vec<_>>();
            let request = assembly_intent::Request {
                material_bindings: Vec::new(),
                name: draft.name.clone(),
                circular: draft.circular,
                fragments: roles
                    .iter()
                    .zip(&draft.slots)
                    .map(
                        |(r, s)| collagene_core::cloning::assembly::intent::Contribution {
                            role: r.role.clone(),
                            input_id: r.input_id.clone(),
                            selector: s.selector.clone(),
                        },
                    )
                    .collect(),
            };
            let intended = if edited {
                assembly_junction_edit::build_intent(&inputs, &request, &draft.junction_edits, None)
                    .unwrap()
                    .target
            } else {
                assembly_intent::build(&inputs, &request, None)
                    .unwrap()
                    .target
            };
            let (catalog, row) = active(&draft).unwrap();
            let native = design::resolve(&catalog, &row, &frame).unwrap();
            assert_eq!(
                native.target.revision(),
                intended.revision(),
                "all name/role/source metadata as well as DNA must match the independently frozen target"
            );
            assert_eq!(native.target, intended);
        }
    }
    #[test]
    fn tailing_rows_keep_explicit_actions_and_scarred_targets_through_projection_and_restart() {
        use collagene_core::{
            cloning::ligation::{
                tailing::TailingBase,
                tailing_recipe::{Action, Profile, ProfileStatus},
            },
            verification::evidence::Status,
        };
        let originals = vec![
            Molecule::dna("Left", "AAACCC", Topology::Linear).unwrap(),
            Molecule::dna("Right", "GGGTTT", Topology::Linear).unwrap(),
        ];
        let frame =
            Workspace::from_documents("Actual tailing sources", originals.clone(), 0).unwrap();
        let mut draft = ordinary();
        draft.circular = true;
        draft.methods.mode = Mode::Tailing;
        for (i, slot) in draft.slots.iter_mut().enumerate() {
            let doc = &frame.documents[i];
            slot.captured = Some(Arc::new(CapturedSource {
                source: SourceInfo {
                    binding: EventBinding {
                        document_id: doc.id.clone(),
                        event_id: doc.head_event.clone().unwrap(),
                        revision: doc.head.clone(),
                    },
                    name: originals[i].name.clone(),
                    length: 6,
                    alphabet: originals[i].sequence.alphabet(),
                    topology: Topology::Linear,
                },
                molecule: originals[i].clone(),
            }));
            slot.selector = collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                start: 0,
                end: 6,
                reverse: false,
            };
        }
        draft.methods.vector.action = Action::AddTail {
            base: TailingBase::T,
        };
        draft.methods.insert.action = Action::AddTail {
            base: TailingBase::A,
        };
        draft.junction_edits = vec![
            Edit {
                junction: 0,
                remove_left: 0,
                remove_right: 0,
                inserted: "T".into(),
            },
            Edit {
                junction: 1,
                remove_left: 0,
                remove_right: 0,
                inserted: "A".into(),
            },
        ];
        draft.sync_rows();
        let first = alternatives::enable(&draft, "ta-library").unwrap();
        let ids = first.alternatives.as_ref().unwrap().rows[0].choices.clone();
        let mut second = alternatives::add_row(&first, ids).unwrap();
        second.methods.tailing_profile = Profile::Gc;
        second.methods.vector.action = Action::AddTail {
            base: TailingBase::C,
        };
        second.methods.insert.action = Action::AddTail {
            base: TailingBase::G,
        };
        second.junction_edits.clear();
        let stored = crate::assembly::draft::model::encode(&second).unwrap();
        let restored = crate::assembly::draft::model::decode(&stored).unwrap();
        let before = crate::assembly::draft::model::encode(&restored).unwrap();
        let projected = project(&restored).unwrap();
        let catalog = projected.request.design.as_ref().unwrap();
        let first = design::resolve(catalog, ready(&projected.request, 0), &frame).unwrap();
        assert_eq!(first.target.sequence.as_str(), "AAACCCTGGGTTTA");
        let second = design::resolve(catalog, ready(&projected.request, 1), &frame).unwrap();
        assert_eq!(
            second.target.sequence.as_str(),
            "AAACCCGGGTTT",
            "missing independently declared scars cannot be supplied by tailing predictions"
        );
        let design::Recipe::Tailing { recipe } = second.recipe else {
            panic!("tailing recipe")
        };
        assert!(recipe.profile == Profile::Gc);
        assert!(matches!(
            recipe.fragments[0].action,
            Action::AddTail {
                base: TailingBase::C
            }
        ));
        let context = batch::begin_revision(&frame, &projected.request).unwrap();
        let mut session = batch::Session::begin(frame.clone(), projected.request.clone()).unwrap();
        let first = session.advance(&context, 0).unwrap().row.unwrap();
        let second = session.advance(&context, 1).unwrap().row.unwrap();
        assert_eq!(first.status, Some(Status::ConsistentUnderAssumptions));
        assert_eq!(second.status, Some(Status::Invalid));
        assert_eq!(
            first.tailing_profile.unwrap().status,
            ProfileStatus::Matched
        );
        assert_eq!(
            second.tailing_profile.unwrap().status,
            ProfileStatus::Matched
        );
        assert_eq!(
            crate::assembly::draft::model::encode(&restored).unwrap(),
            before
        );
    }
}

#[cfg(test)]
#[path = "tests/paired_oligo_library.rs"]
mod paired_oligo_library;
