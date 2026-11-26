//! Shared conversion of editable preparation rows to existing typed recipes.
use crate::cloning_methods_model::{MethodsState, OverlapPreparation};
use collagene_core::{Molecule, analysis::Primer, verification::PlanInput};

pub(crate) fn prepared_rows(
    state: &MethodsState,
    input_count: usize,
) -> Result<
    (
        Vec<String>,
        Vec<collagene_core::assembly_recipe::FragmentSpec>,
    ),
    String,
> {
    use collagene_core::assembly_recipe::{FragmentSpec, Preparation, RunPreparation};
    if state.rows.is_empty() || state.rows.len() > 12 || input_count > 24 {
        return Err("Choose 1–12 ordered preparation rows and their exact original sources".into());
    }
    let roles = state.role_ids(state.rows.len())?;
    let ids = state.all_input_ids()?;
    if ids.len() != input_count {
        return Err("Preparation inputs differ from the complete pinned source list.".into());
    }
    let mut fragments = Vec::with_capacity(state.rows.len());
    for (i, row) in state.rows.iter().enumerate() {
        if let Some(request) = &row.overlap_extension_pcr {
            request.validate_raw()?;
        }
        if row.fragment.len() > 20 || row.forward_name.len() > 256 || row.reverse_name.len() > 256 {
            return Err("Preparation names or fragment number exceed their field bounds".into());
        }
        if [
            &row.enzymes,
            &row.forward,
            &row.reverse_primer,
            &row.forward_tail,
            &row.reverse_tail,
        ]
        .iter()
        .any(|s| s.len() > 65536)
        {
            return Err("Preparation fields exceed the 64 KiB field limit".into());
        }
        let preparation = match row.preparation {
            OverlapPreparation::AnnealedOligos=>Preparation::AnnealedOligos {
                bottom_input_id:state.bottom_input_id(i)?,
                options:row.anneal_options.clone(),
                candidate_index:row.anneal_candidate.ok_or("Inspect the annealing alternatives and explicitly choose a candidate before assembling.")?,
            },
            OverlapPreparation::WholeLinear => Preparation::WholeLinear,
            OverlapPreparation::OverlapExtensionPcr => Preparation::OverlapExtensionPcr {
                request: row.overlap_extension_pcr.as_ref().ok_or("Enter the two source-PCR pairs, terminal homology and outer-PCR pair before reconstruction.")?.request()?,
            },
            OverlapPreparation::OverlapExtensionPcrCatalogDigest => {
                let request = row.overlap_extension_pcr.as_ref().ok_or("Enter all six overlap-extension oligos before inspecting the outer-PCR digest.")?.request()?;
                let selection = row.catalog_selection.as_ref().ok_or("Pin the outer-PCR product digest enzymes first.")?;
                selection.validate().map_err(|e| e.to_string())?;
                let fragment_index = row.fragment.trim().parse::<usize>().ok().and_then(|n| n.checked_sub(1)).ok_or("Digest fragment numbers start at1")?;
                Preparation::OverlapExtensionPcrCatalogDigest { request, selection: selection.clone(), fragment_index, repair: row.repair }
            }
            OverlapPreparation::DeclaredDuplex => Preparation::DeclaredDuplex {
                left: row.left_terminus,
                right: row.right_terminus,
            },
            OverlapPreparation::CatalogDigest | OverlapPreparation::PcrCatalogDigest | OverlapPreparation::MismatchPcrCatalogDigest => {
                let selection = row.catalog_selection.as_ref().ok_or(
                    "Pin a catalog selection using Virtual digest or the PCR-product enzyme selector; then choose its fragment.",
                )?;
                selection.validate().map_err(|e| e.to_string())?;
                let fragment_index = row
                    .fragment
                    .trim()
                    .parse::<usize>()
                    .ok()
                    .and_then(|n| n.checked_sub(1))
                    .ok_or("Digest fragment numbers start at 1")?;
                if row.preparation == OverlapPreparation::MismatchPcrCatalogDigest {
                    Preparation::MismatchPcrCatalogDigest {
                        forward: Primer::new(&row.forward_name, &row.forward, &row.forward_tail).map_err(|e|e.to_string())?,
                        reverse: Primer::new(&row.reverse_name, &row.reverse_primer, &row.reverse_tail).map_err(|e|e.to_string())?,
                        options: row.mismatch_options()?, selection: selection.clone(), fragment_index, repair: row.repair,
                    }
                } else if row.preparation == OverlapPreparation::PcrCatalogDigest {
                    Preparation::PcrCatalogDigest {
                        forward: Primer::new(&row.forward_name, &row.forward, &row.forward_tail)
                            .map_err(|e| e.to_string())?,
                        reverse: Primer::new(
                            &row.reverse_name,
                            &row.reverse_primer,
                            &row.reverse_tail,
                        )
                        .map_err(|e| e.to_string())?,
                        selection: selection.clone(),
                        fragment_index,
                        repair: row.repair,
                    }
                } else {
                    Preparation::CatalogDigest {
                        selection: selection.clone(),
                        fragment_index,
                        repair: row.repair,
                    }
                }
            }
            OverlapPreparation::Digest => {
                let enzymes = row
                    .enzymes
                    .split(',')
                    .map(str::trim)
                    .filter(|s| !s.is_empty())
                    .map(str::to_owned)
                    .collect::<Vec<_>>();
                if enzymes.is_empty() || enzymes.len() > 64 || enzymes.iter().any(|e| e.len() > 256)
                {
                    return Err("Choose 1–64 explicit restriction enzyme names".into());
                }
                let fragment_index = row
                    .fragment
                    .trim()
                    .parse::<usize>()
                    .ok()
                    .and_then(|n| n.checked_sub(1))
                    .ok_or("Digest fragment numbers start at 1")?;
                Preparation::Digest {
                    enzymes,
                    fragment_index,
                    repair: row.repair,
                }
            }
            OverlapPreparation::MismatchPcr => Preparation::MismatchPcr {
                forward: Primer::new(&row.forward_name, &row.forward, &row.forward_tail).map_err(|e|e.to_string())?,
                reverse: Primer::new(&row.reverse_name, &row.reverse_primer, &row.reverse_tail).map_err(|e|e.to_string())?,
                options: row.mismatch_options()?,
            },
            OverlapPreparation::ExactPcr => Preparation::ExactPcr {
                forward: Primer::new(&row.forward_name, &row.forward, &row.forward_tail)
                    .map_err(|e| e.to_string())?,
                reverse: Primer::new(&row.reverse_name, &row.reverse_primer, &row.reverse_tail)
                    .map_err(|e| e.to_string())?,
            },
        };
        let preparation = if let Some(origin) = &row.run_origin {
            origin.validate_documentary().map_err(|e| e.to_string())?;
            let preparation = match preparation {
                Preparation::WholeLinear => RunPreparation::WholeLinear,
                Preparation::Digest { enzymes, fragment_index, repair } => RunPreparation::Digest { enzymes, fragment_index, repair },
                Preparation::CatalogDigest { selection, fragment_index, repair } => RunPreparation::CatalogDigest { selection, fragment_index, repair },
                Preparation::ExactPcr { forward, reverse } => RunPreparation::ExactPcr { forward, reverse },
                Preparation::OverlapExtensionPcr { request } => RunPreparation::OverlapExtensionPcr { request },
                Preparation::OverlapExtensionPcrCatalogDigest { request, selection, fragment_index, repair } => RunPreparation::OverlapExtensionPcrCatalogDigest { request, selection, fragment_index, repair },
                Preparation::MismatchPcr { forward, reverse, options } => RunPreparation::MismatchPcr { forward, reverse, options },
                Preparation::PcrCatalogDigest { forward, reverse, selection, fragment_index, repair } => RunPreparation::PcrCatalogDigest { forward, reverse, selection, fragment_index, repair },
                Preparation::MismatchPcrCatalogDigest { forward, reverse, options, selection, fragment_index, repair } => RunPreparation::MismatchPcrCatalogDigest { forward, reverse, options, selection, fragment_index, repair },
                Preparation::DeclaredDuplex { .. } | Preparation::AnnealedOligos { .. } | Preparation::FromRun { .. } => return Err("Saved-run material supports whole material, digest or exact PCR preparation; declared duplex and paired-oligo substitution are unsupported.".into()),
            };
            Preparation::FromRun {
                origin: origin.as_ref().clone(),
                preparation,
            }
        } else {
            preparation
        };
        fragments.push(FragmentSpec {
            role: roles[i].role.clone(),
            input_id: ids[i].clone(),
            preparation,
            reverse: row.reverse,
        });
    }
    Ok((ids, fragments))
}

pub(crate) fn tailing_request(
    state: &MethodsState,
    inputs: Vec<Molecule>,
    target: &Molecule,
) -> Result<(Vec<PlanInput>, collagene_core::tailing_recipe::Recipe), String> {
    use collagene_core::tailing_recipe::{Fragment, Recipe};
    if inputs.len() != 2 {
        return Err("TA/GC cloning requires exactly two ordered source fragments.".into());
    }
    let (ids, sources) = prepared_rows(state, inputs.len())?;
    let recipe = Recipe {
        name: format!("Declared tailing and ligation to {}", target.name),
        fragments: sources
            .into_iter()
            .zip([&state.vector, &state.insert])
            .map(|(source, settings)| Fragment {
                source,
                base: settings.base,
                treatment: settings.treatment.operation(),
            })
            .collect(),
        circular: state.circular,
        allow_nicked: state.allow_nicked,
    };
    let borrowed = ids
        .iter()
        .zip(&inputs)
        .map(|(id, m)| (id.as_str(), m))
        .collect::<Vec<_>>();
    collagene_core::tailing_recipe::source_preflight(&borrowed, target, &recipe)
        .map_err(|e| e.to_string())?;
    let bound = ids
        .into_iter()
        .zip(inputs)
        .map(|(id, molecule)| PlanInput {
            revision: molecule.revision(),
            id,
            molecule,
        })
        .collect();
    Ok((bound, recipe))
}

#[cfg(test)]
mod tailing_tests {
    use super::*;
    use crate::{
        assembly_role_ids::Role,
        cloning_methods_model::{Mode, Treatment},
    };
    use collagene_core::{Topology, chemistry::EndTreatment, tailing::TailingBase};
    #[test]
    fn pinned_tailing_request_uses_rows_and_exact_role_ids() {
        let mut state = MethodsState::default();
        state.mode = Mode::Tailing;
        state.source_roles = Some(vec![
            Role {
                role: "backbone".into(),
                input_id: "captured-backbone".into(),
            },
            Role {
                role: "cargo".into(),
                input_id: "captured-cargo".into(),
            },
        ]);
        state.rows[0].preparation = OverlapPreparation::WholeLinear;
        state.rows[1].preparation = OverlapPreparation::WholeLinear;
        state.rows[1].reverse = true;
        state.vector.base = TailingBase::T;
        state.insert.base = TailingBase::A;
        state.insert.treatment = Treatment::Phosphorylate;
        // Legacy digest text must not replace row preparations or captured IDs.
        state.vector.enzymes = "Legacy setting is deliberately unusable".into();
        let originals = vec![
            Molecule::dna("backbone", "AACC", Topology::Linear).unwrap(),
            Molecule::dna("cargo", "GGTT", Topology::Linear).unwrap(),
        ];
        let target = Molecule::dna("independent intent", "AACCTAACCA", Topology::Circular).unwrap();
        let (bound, recipe) = tailing_request(&state, originals.clone(), &target).unwrap();
        assert_eq!(bound[0].id, "captured-backbone");
        assert_eq!(bound[1].id, "captured-cargo");
        assert_eq!(bound[0].revision, originals[0].revision());
        assert_eq!(bound[1].molecule, originals[1]);
        assert_eq!(recipe.fragments[1].source.role, "cargo");
        assert!(recipe.fragments[1].source.reverse);
        assert!(matches!(
            recipe.fragments[0].source.preparation,
            collagene_core::assembly_recipe::Preparation::WholeLinear
        ));
        assert_eq!(recipe.fragments[0].base, TailingBase::T);
        assert_eq!(
            recipe.fragments[1].treatment,
            Some(EndTreatment::Phosphorylate)
        );
        assert_eq!(target.sequence.as_str(), "AACCTAACCA");
    }
    #[test]
    fn independently_declared_circle_and_linear_terminal_intent_feed_actual_tailing() {
        use collagene_core::{
            assembly_intent, assembly_junction_edit, overlap_primer_design, verification::Status,
        };
        for circular in [true, false] {
            let mut state = MethodsState::default();
            state.mode = Mode::Tailing;
            state.circular = circular;
            for (i, (forward, reverse)) in [("AAA", "GGG"), ("GGG", "AAA")].into_iter().enumerate()
            {
                state.rows[i].preparation = OverlapPreparation::ExactPcr;
                state.rows[i].forward = forward.into();
                state.rows[i].reverse_primer = reverse.into();
            }
            state.vector.base = TailingBase::T;
            state.insert.base = TailingBase::A;
            state.vector.treatment = Treatment::Phosphorylate;
            state.insert.treatment = Treatment::Phosphorylate;
            let originals = vec![
                Molecule::dna("backbone", "AAACCC", Topology::Linear).unwrap(),
                Molecule::dna("cargo", "GGGTTT", Topology::Linear).unwrap(),
            ];
            let roles = state.role_ids(2).unwrap();
            let inputs = originals
                .iter()
                .zip(&roles)
                .map(|(m, role)| PlanInput {
                    id: role.input_id.clone(),
                    revision: m.revision(),
                    molecule: m.clone(),
                })
                .collect::<Vec<_>>();
            let base = assembly_intent::Request {
                name: "Operator declared target".into(),
                circular,
                fragments: roles
                    .iter()
                    .map(|role| overlap_primer_design::Fragment {
                        role: role.role.clone(),
                        input_id: role.input_id.clone(),
                        start: 0,
                        end: 6,
                        reverse: false,
                    })
                    .collect(),
            };
            let edits = vec![
                assembly_junction_edit::Edit {
                    junction: 0,
                    remove_left: 0,
                    remove_right: 0,
                    inserted: "T".into(),
                },
                assembly_junction_edit::Edit {
                    junction: 1,
                    remove_left: 0,
                    remove_right: 0,
                    inserted: "A".into(),
                },
            ];
            let intent = assembly_junction_edit::build_intent(&inputs, &base, &edits).unwrap();
            assert_eq!(intent.target.sequence.as_str(), "AAACCCTGGGTTTA");
            let (bound, recipe) = tailing_request(&state, originals, &intent.target).unwrap();
            let preview =
                collagene_core::tailing_recipe::preview(&bound, &intent.target, &recipe).unwrap();
            assert_eq!(
                preview.verification.status,
                Status::ConsistentUnderAssumptions,
                "{:?}",
                preview.verification.findings
            );
            assert_eq!(
                preview.verification.products.last().unwrap().sequence.len(),
                14
            );
            assert_eq!(
                preview.verification.products.last().unwrap().topology,
                if circular {
                    Topology::Circular
                } else {
                    Topology::Linear
                }
            );
            assert_eq!(preview.junctions.len(), if circular { 2 } else { 1 });
            assert_eq!(intent.target.sequence.as_str(), "AAACCCTGGGTTTA");
        }
    }
    #[test]
    fn tailing_request_rejects_missing_or_extra_sources() {
        let state = MethodsState::default();
        let target = Molecule::dna("intent", "AACC", Topology::Linear).unwrap();
        assert!(tailing_request(&state, vec![], &target).is_err());
        assert!(tailing_request(&state, vec![target.clone(); 3], &target).is_err());
    }
}

#[cfg(test)]
mod annealed_tests {
    use super::*;
    use crate::cloning_methods_model::{Mode, Treatment};
    use collagene_core::{
        Topology, assembly_intent, catalog_digest, overlap_primer_design, verification::Status,
    };
    #[test]
    fn two_reference_contributions_keep_three_original_sources_through_ligation() {
        let catalog = catalog_digest::builtin_catalog().unwrap();
        let ecori = catalog
            .edition()
            .enzymes
            .iter()
            .find(|e| e.name == "EcoRI")
            .unwrap()
            .id
            .clone();
        let mut state = MethodsState::default();
        state.mode = Mode::Ligation;
        state.circular = true;
        state.rows[0].preparation = OverlapPreparation::CatalogDigest;
        state.rows[0].catalog_selection = Some(catalog_digest::select(&catalog, &[ecori]).unwrap());
        state.rows[1].preparation = OverlapPreparation::AnnealedOligos;
        state.rows[1].anneal_bottom_document = Some(2);
        state.rows[1].anneal_options.min_paired_bases = 6;
        state.rows[1].anneal_candidate = Some(0);
        state.rows[1].treatment = Treatment::Phosphorylate;
        let originals = vec![
            Molecule::dna("Backbone", "GGAATTCTTACG", Topology::Circular).unwrap(),
            Molecule::dna("Top", "AATTGCTACG", Topology::Linear).unwrap(),
            Molecule::dna("Bottom", "AATTCGTAGC", Topology::Linear).unwrap(),
        ];
        let ids = state.all_input_ids().unwrap();
        let inputs = ids
            .iter()
            .zip(&originals)
            .map(|(id, m)| PlanInput {
                id: id.clone(),
                molecule: m.clone(),
                revision: m.revision(),
            })
            .collect::<Vec<_>>();
        let roles = state.role_ids(2).unwrap();
        let declaration = assembly_intent::Request {
            name: "Independent backbone and top contribution".into(),
            circular: true,
            fragments: vec![
                overlap_primer_design::Fragment {
                    role: roles[0].role.clone(),
                    input_id: ids[0].clone(),
                    start: 2,
                    end: 14,
                    reverse: false,
                },
                overlap_primer_design::Fragment {
                    role: roles[1].role.clone(),
                    input_id: ids[1].clone(),
                    start: 0,
                    end: 10,
                    reverse: false,
                },
            ],
        };
        let intent = assembly_intent::build(&inputs, &declaration).unwrap();
        assert_eq!(intent.target.sequence.as_str(), "AATTCTTACGGGAATTGCTACG");
        let (bound_ids, fragments) = prepared_rows(&state, 3).unwrap();
        assert_eq!(bound_ids, ids);
        assert_eq!(fragments.len(), 2);
        assert!(
            matches!(&fragments[1].preparation,collagene_core::assembly_recipe::Preparation::AnnealedOligos{bottom_input_id,..} if bottom_input_id==&ids[2])
        );
        let recipe = collagene_core::ligation_recipe::Recipe {
            name: "Actual three-source ligation".into(),
            fragments: fragments
                .into_iter()
                .zip(&state.rows)
                .map(|(source, row)| collagene_core::ligation_recipe::Fragment {
                    source,
                    treatment: row.treatment.operation(),
                })
                .collect(),
            circular: true,
            allow_nicked: false,
        };
        let preview =
            collagene_core::ligation_recipe::preview(&inputs, &intent.target, &recipe).unwrap();
        assert_eq!(preview.compiled.plan.inputs.len(), 3);
        assert_eq!(preview.compiled.rows.len(), 2);
        assert_eq!(
            preview.verification.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            preview.verification.findings
        );
        assert_eq!(
            preview
                .verification
                .products
                .last()
                .unwrap()
                .sequence
                .as_str(),
            intent.target.sequence.as_str()
        );
        state.rows[1].anneal_candidate = None;
        assert!(
            prepared_rows(&state, 3)
                .unwrap_err()
                .contains("explicitly choose")
        );
        assert!(prepared_rows(&state, 2).is_err());
    }
}

#[cfg(test)]
mod mismatch_tests {
    use super::*;
    #[test]
    fn mismatch_preparation_retains_named_oligo_domains_policy_and_source_role() {
        let mut state = MethodsState::default();
        state.rows.truncate(1);
        state.source_roles = Some(vec![crate::assembly_role_ids::Role {
            role: "cargo".into(),
            input_id: "exact-cargo-event".into(),
        }]);
        let row = &mut state.rows[0];
        row.preparation = OverlapPreparation::MismatchPcr;
        row.forward_name = "Cargo mutation F".into();
        row.reverse_name = "Cargo mutation R".into();
        row.forward = "ACGACCCC".into();
        row.reverse_primer = "GTTACCCC".into();
        row.forward_tail = "GG".into();
        row.reverse_tail = "TT".into();
        row.mismatch_max_mismatches = "2".into();
        row.mismatch_three_prime_exact = "3".into();
        let before = state.clone();
        let (ids, rows) = prepared_rows(&state, 1).unwrap();
        assert_eq!(ids, vec!["exact-cargo-event"]);
        assert_eq!(rows[0].role, "cargo");
        let collagene_core::assembly_recipe::Preparation::MismatchPcr {
            forward,
            reverse,
            options,
        } = &rows[0].preparation
        else {
            panic!("explicit mismatch preparation")
        };
        assert_eq!(forward.name, "Cargo mutation F");
        assert_eq!(forward.ordered_sequence(), "GGACGACCCC");
        assert_eq!(reverse.ordered_sequence(), "TTGTTACCCC");
        assert_eq!(options.binding.max_mismatches, 2);
        assert_eq!(options.binding.min_three_prime_exact, 3);
        assert_eq!(
            options.max_pair_checks,
            collagene_core::mutagenesis::PcrOptions::default().max_pair_checks
        );
        assert!(state == before);
        state.rows[0].mismatch_max_mismatches = "99".into();
        assert!(prepared_rows(&state, 1).is_err());
    }
}

#[cfg(test)]
mod overlap_extension_tests {
    use super::*;
    use crate::cloning_methods_model::{ExtensionPrimerInput, Mode, OverlapExtensionPcrInput};
    use collagene_core::{
        Topology, assembly_recipe, cloning_outputs,
        run_lineage::{EventBinding, InputEventBinding},
        run_material,
        verification::{DesignPlan, Operation, Source, Status},
        workspace::Workspace,
    };
    const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
    const EDITED: &str = "AGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCA";
    const FINAL: &str = "AGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCAGGCTAACCGTTA";
    fn primer(name: &str, annealing: &str, tail: &str) -> ExtensionPrimerInput {
        ExtensionPrimerInput {
            name: name.into(),
            annealing: annealing.into(),
            tail: tail.into(),
        }
    }
    fn raw() -> OverlapExtensionPcrInput {
        let mut r = OverlapExtensionPcrInput::default();
        r.left.forward = primer("Left outer", "AGTCCGAT", "");
        r.left.reverse = primer("Left mutagenic", "TACGTAGC", "GGTCTGAAAA");
        r.right.forward = primer("Right mutagenic", "TTCAGACC", "GCTACGTATT");
        r.right.reverse = primer("Right outer", "TGACGATC", "");
        r.outer.forward = primer("Final forward", "AGTCCGAT", "");
        r.outer.reverse = primer("Final reverse", "TGACGATC", "");
        r.overlap = "18".into();
        r.min_paired_bases = "8".into();
        r
    }
    fn state() -> MethodsState {
        let mut s = MethodsState::default();
        s.mode = Mode::Overlap;
        s.circular = false;
        s.overlaps = "8".into();
        s.source_roles = Some(vec![
            crate::assembly_role_ids::Role {
                role: "edited-insert".into(),
                input_id: "original-stock-event".into(),
            },
            crate::assembly_role_ids::Role {
                role: "second-part".into(),
                input_id: "other-event".into(),
            },
        ]);
        s.rows[0].preparation = OverlapPreparation::OverlapExtensionPcr;
        s.rows[0].overlap_extension_pcr = Some(raw());
        // The previous single-pair fields are deliberately unusable, and must not
        // replace any of the six new oligos during conversion.
        s.rows[0].forward = "wrong old primer".into();
        s.rows[0].reverse_primer = "wrong old primer".into();
        s
    }
    fn preview(s: &MethodsState, source: Molecule) -> assembly_recipe::Preview {
        let (ids, fragments) = prepared_rows(s, 2).unwrap();
        let target =
            Molecule::dna("Independent two-part requirement", FINAL, Topology::Linear).unwrap();
        let other = Molecule::dna(
            "Unchanged second part",
            "GATCGTCAGGCTAACCGTTA",
            Topology::Linear,
        )
        .unwrap();
        let inputs = ids
            .into_iter()
            .zip([source, other])
            .map(|(id, molecule)| PlanInput {
                id,
                revision: molecule.revision(),
                molecule,
            })
            .collect::<Vec<_>>();
        let (_, options) = s.profile().unwrap();
        let recipe = assembly_recipe::Recipe {
            name: "Actual edited insert and separate part".into(),
            fragments,
            overlaps: vec![8],
            circular: false,
            options,
        };
        assembly_recipe::preview(&inputs, &target, &recipe).unwrap()
    }
    #[test]
    fn six_manual_oligos_map_into_real_multipart_preparation_without_changing_original() {
        let s = state();
        let before = s.clone();
        let original =
            Molecule::dna("Original unchanged stock", ORIGINAL, Topology::Linear).unwrap();
        let original_json = serde_json::to_value(&original).unwrap();
        let p = preview(&s, original.clone());
        assert_eq!(
            p.verification.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            p.verification.findings
        );
        assert_eq!(p.compiled.plan.inputs[0].id, "original-stock-event");
        assert_eq!(p.compiled.rows[0].role, "edited-insert");
        assert_eq!(p.prepared[0].sequence.as_str(), EDITED);
        assert_eq!(
            p.verification.products[0].sequence.as_str(),
            "AGTCCGATGCTACGTATTTTCAGACC"
        );
        assert_eq!(
            p.verification.products[1].sequence.as_str(),
            "GCTACGTATTTTCAGACCTGTAGATCGTCA"
        );
        assert_eq!(p.verification.products[2].sequence.as_str(), EDITED);
        assert_eq!(p.verification.products[3].sequence.as_str(), EDITED);
        assert_eq!(
            p.verification.products.last().unwrap().sequence.as_str(),
            FINAL
        );
        assert_eq!(
            serde_json::to_value(&p.compiled.plan.inputs[0].molecule).unwrap(),
            original_json
        );
        assert!(s == before);
        let mut incomplete = s;
        incomplete.rows[0]
            .overlap_extension_pcr
            .as_mut()
            .unwrap()
            .outer
            .reverse
            .annealing
            .clear();
        assert!(prepared_rows(&incomplete, 2).is_err());
    }
    #[test]
    fn saved_run_child_keeps_origin_and_uses_both_pcr_branches_on_reconstructed_material() {
        let original = Molecule::dna("Original", ORIGINAL, Topology::Linear).unwrap();
        let target = Molecule::dna(
            "Separate original PCR requirement",
            ORIGINAL,
            Topology::Linear,
        )
        .unwrap();
        let w = Workspace::from_documents(
            "Original preparation",
            vec![original.clone(), target.clone()],
            0,
        )
        .unwrap();
        let event = |i: usize| {
            let d = &w.documents[i];
            EventBinding {
                document_id: d.id.clone(),
                event_id: d.head_event.clone().unwrap(),
                revision: d.head.clone(),
            }
        };
        let plan = DesignPlan {
            name: "Actual original-stock PCR".into(),
            target_revision: target.revision(),
            inputs: vec![PlanInput {
                id: "stock".into(),
                revision: original.revision(),
                molecule: original,
            }],
            steps: vec![Operation::Pcr {
                template: Source::Input("stock".into()),
                forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
                reverse: Primer::new("R", "TGACGATC", "").unwrap(),
            }],
        };
        let commit = cloning_outputs::finalize(
            &w,
            &plan,
            event(1),
            vec![InputEventBinding {
                input_id: "stock".into(),
                event: event(0),
            }],
            &cloning_outputs::Options {
                name: "Saved actual stock PCR".into(),
                save_primers: false,
                show_prepared: false,
                allow_incomplete: false,
            },
        )
        .unwrap();
        let selection = run_material::Selection {
            run_id: commit.receipt.run.run_id.clone(),
            step_index: 0,
            event: commit.receipt.construct.clone(),
        };
        let capsule = run_material::capsule(&commit.workspace, &selection).unwrap();
        let saved = commit.workspace.revisions[&selection.event.revision].clone();
        let mut s = state();
        s.rows[0].run_origin =
            Some(crate::run_origin_model::PinnedOrigin::new(capsule.clone()).unwrap());
        let (_, rows) = prepared_rows(&s, 2).unwrap();
        let assembly_recipe::Preparation::FromRun {
            origin,
            preparation: assembly_recipe::RunPreparation::OverlapExtensionPcr { request },
        } = &rows[0].preparation
        else {
            panic!("flat original-run child");
        };
        assert_eq!(
            serde_json::to_value(origin).unwrap(),
            serde_json::to_value(&capsule).unwrap()
        );
        assert_eq!(
            request.left.reverse.ordered_sequence(),
            "GGTCTGAAAATACGTAGC"
        );
        assert_eq!(
            request.right.forward.ordered_sequence(),
            "GCTACGTATTTTCAGACC"
        );
        let p = preview(&s, saved);
        assert_eq!(
            p.verification.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            p.verification.findings
        );
        assert!(matches!(
            &p.compiled.plan.steps[0],
            Operation::ReconstructRunProduct { .. }
        ));
        assert!(matches!(
            &p.compiled.plan.steps[1],
            Operation::Pcr {
                template: Source::Step(0),
                ..
            }
        ));
        assert!(matches!(
            &p.compiled.plan.steps[2],
            Operation::Pcr {
                template: Source::Step(0),
                ..
            }
        ));
        assert_eq!(
            p.verification.products.last().unwrap().sequence.as_str(),
            FINAL
        );
        s.rows[0].preparation = OverlapPreparation::OverlapExtensionPcrCatalogDigest;
        let catalog = collagene_core::catalog_digest::builtin_catalog().unwrap();
        let id = catalog
            .edition()
            .enzymes
            .iter()
            .find(|entry| entry.name == "EcoRI")
            .unwrap()
            .id
            .clone();
        s.rows[0].catalog_selection =
            Some(collagene_core::catalog_digest::select(&catalog, &[id]).unwrap());
        s.rows[0].fragment = "2".into();
        let (_, rows) = prepared_rows(&s, 2).unwrap();
        let assembly_recipe::Preparation::FromRun {
            origin,
            preparation:
                assembly_recipe::RunPreparation::OverlapExtensionPcrCatalogDigest {
                    request,
                    selection,
                    fragment_index,
                    ..
                },
        } = &rows[0].preparation
        else {
            panic!("saved-source OE digest child");
        };
        assert_eq!(
            serde_json::to_value(origin).unwrap(),
            serde_json::to_value(&capsule).unwrap()
        );
        assert_eq!(
            request.left.reverse.ordered_sequence(),
            "GGTCTGAAAATACGTAGC"
        );
        assert_eq!(selection, s.rows[0].catalog_selection.as_ref().unwrap());
        assert_eq!(*fragment_index, 1);
    }
    #[test]
    fn extension_catalog_conversion_preserves_six_oligos_exact_pin_and_repair_request() {
        let mut s = state();
        let row = &mut s.rows[0];
        row.preparation = OverlapPreparation::OverlapExtensionPcrCatalogDigest;
        row.overlap_extension_pcr
            .as_mut()
            .unwrap()
            .outer
            .forward
            .tail = "CCGAATTC".into();
        row.overlap_extension_pcr
            .as_mut()
            .unwrap()
            .outer
            .reverse
            .tail = "GGCTCGAG".into();
        row.fragment = "2".into();
        row.repair = Some(collagene_core::chemistry::EndRepairMethod::Blunt);
        row.reverse = true;
        let catalog = collagene_core::catalog_digest::builtin_catalog().unwrap();
        let ids = catalog
            .edition()
            .enzymes
            .iter()
            .filter(|entry| entry.name == "EcoRI" || entry.name == "XhoI")
            .map(|entry| entry.id.clone())
            .collect::<Vec<_>>();
        row.catalog_selection =
            Some(collagene_core::catalog_digest::select(&catalog, &ids).unwrap());
        let before = s.clone();
        let (_, rows) = prepared_rows(&s, 2).unwrap();
        let assembly_recipe::Preparation::OverlapExtensionPcrCatalogDigest {
            request,
            selection,
            fragment_index,
            repair,
        } = &rows[0].preparation
        else {
            panic!("OE outer-PCR catalog digest");
        };
        assert_eq!(
            request.left.reverse.ordered_sequence(),
            "GGTCTGAAAATACGTAGC"
        );
        assert_eq!(
            request.right.forward.ordered_sequence(),
            "GCTACGTATTTTCAGACC"
        );
        assert_eq!(request.outer.forward.tail, "CCGAATTC");
        assert_eq!(request.outer.reverse.tail, "GGCTCGAG");
        assert_eq!(selection, s.rows[0].catalog_selection.as_ref().unwrap());
        assert_eq!(*fragment_index, 1);
        assert_eq!(
            *repair,
            Some(collagene_core::chemistry::EndRepairMethod::Blunt)
        );
        assert!(rows[0].reverse);
        assert!(s == before);
        s.rows[0].catalog_selection = None;
        assert!(prepared_rows(&s, 2).is_err());
    }
}
