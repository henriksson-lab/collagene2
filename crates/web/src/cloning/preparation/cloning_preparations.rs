//! Shared conversion of editable preparation rows to existing typed recipes.
use crate::cloning::preparation::cloning_methods_model::OverlapPreparation;
use crate::cloning::preparation::cloning_methods_model::{MethodsState, Mode};
use collagene_core::{planning::plan::PlanInput, sequence::analysis::Primer};
use collagene_model::Molecule;

pub(crate) fn prepared_rows(
    state: &MethodsState,
    input_count: usize,
) -> Result<
    (
        Vec<String>,
        Vec<collagene_core::cloning::assembly::recipe::FragmentSpec>,
    ),
    String,
> {
    use collagene_core::cloning::assembly::recipe::{FragmentSpec, Preparation, RunPreparation};
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
        row.reaction_conditions.validate_raw()?;
        let conditions = if row.reaction_conditions.enabled {
            if !matches!(
                row.preparation,
                OverlapPreparation::CatalogDigest
                    | OverlapPreparation::PcrCatalogDigest
                    | OverlapPreparation::MismatchPcrCatalogDigest
                    | OverlapPreparation::OverlapExtensionPcrCatalogDigest
            ) {
                return Err("Reaction conditions require a catalog-pinned digest preparation. Disable the retained setup or choose that preparation.".into());
            }
            let selection = row
                .catalog_selection
                .as_ref()
                .ok_or("Pin the exact digest enzymes before enabling reaction conditions.")?;
            row.reaction_conditions.options(
                &collagene_core::cloning::restriction::catalog::builtin_catalog()
                    .map_err(|e| e.to_string())?,
                selection,
            )?
        } else {
            None
        };
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
                Preparation::OverlapExtensionPcrCatalogDigest { request, selection: selection.clone(), conditions, fragment_index, repair: row.repair }
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
                        options: row.mismatch_options()?, selection: selection.clone(), conditions, fragment_index, repair: row.repair,
                    }
                } else if row.preparation == OverlapPreparation::PcrCatalogDigest {
                    Preparation::PcrCatalogDigest {
                        template_strand:row.template_strand,
                        forward: Primer::new(&row.forward_name, &row.forward, &row.forward_tail)
                            .map_err(|e| e.to_string())?,
                        reverse: Primer::new(
                            &row.reverse_name,
                            &row.reverse_primer,
                            &row.reverse_tail,
                        )
                        .map_err(|e| e.to_string())?,
                        selection: selection.clone(),
                        conditions,
                        fragment_index,
                        repair: row.repair,
                    }
                } else {
                    Preparation::CatalogDigest {
                        selection: selection.clone(),
                        conditions,
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
                template_strand:row.template_strand,
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
                Preparation::CatalogDigest {
                    conditions, selection, fragment_index, repair } => RunPreparation::CatalogDigest {
                    conditions, selection, fragment_index, repair },
                Preparation::ExactPcr { template_strand,forward, reverse } => RunPreparation::ExactPcr { template_strand,forward, reverse },
                Preparation::OverlapExtensionPcr { request } => RunPreparation::OverlapExtensionPcr { request },
                Preparation::OverlapExtensionPcrCatalogDigest {
                    conditions, request, selection, fragment_index, repair } => RunPreparation::OverlapExtensionPcrCatalogDigest {
                    conditions, request, selection, fragment_index, repair },
                Preparation::MismatchPcr { forward, reverse, options } => RunPreparation::MismatchPcr { forward, reverse, options },
                Preparation::PcrCatalogDigest {
                    conditions, template_strand,forward, reverse, selection, fragment_index, repair } => RunPreparation::PcrCatalogDigest {
                    conditions, template_strand,forward, reverse, selection, fragment_index, repair },
                Preparation::MismatchPcrCatalogDigest {
                    conditions, forward, reverse, options, selection, fragment_index, repair } => RunPreparation::MismatchPcrCatalogDigest {
                    conditions, forward, reverse, options, selection, fragment_index, repair },
                Preparation::DeclaredDuplex { .. } | Preparation::AnnealedOligos { .. } | Preparation::FromRun { .. } | Preparation::FromDerived { .. } => return Err("Saved-run material supports whole material, digest or exact PCR preparation; declared duplex and paired-oligo substitution are unsupported.".into()),
            };
            match origin.as_ref() {
                collagene_core::material::run::MaterialOrigin::Run(origin) => {
                    Preparation::FromRun {
                        origin: origin.clone(),
                        preparation,
                    }
                }
                collagene_core::material::run::MaterialOrigin::Derived(origin) => {
                    Preparation::FromDerived {
                        origin: origin.clone(),
                        preparation,
                    }
                }
            }
        } else {
            preparation
        };
        fragments.push(FragmentSpec {
            generation: row.generation()?,
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
) -> Result<
    (
        Vec<PlanInput>,
        collagene_core::cloning::ligation::tailing_recipe::Recipe,
    ),
    String,
> {
    use collagene_core::cloning::ligation::tailing_recipe::{Fragment, Recipe};
    if state.rows.len() != 2 {
        return Err("TA/GC cloning requires exactly two ordered prepared fragments.".into());
    }
    let (ids, sources) = prepared_rows(state, inputs.len())?;
    let recipe = Recipe {
        name: format!("Declared tailing and ligation to {}", target.name),
        profile: state.tailing_profile,
        fragments: sources
            .into_iter()
            .zip([&state.vector, &state.insert])
            .map(|(source, settings)| Fragment {
                source,
                action: settings.action,
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
    collagene_core::cloning::ligation::tailing_recipe::source_preflight(&borrowed, target, &recipe)
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
    use crate::assembly::draft::assembly_role_ids::Role;
    use crate::cloning::preparation::cloning_methods_model::Mode;
    use crate::cloning::preparation::cloning_methods_model::Treatment;
    use collagene_core::{chemistry::EndTreatment, cloning::ligation::tailing::TailingBase};
    use collagene_model::Topology;
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
        state.vector.action = collagene_core::cloning::ligation::tailing_recipe::Action::AddTail {
            base: TailingBase::T,
        };
        state.insert.action = collagene_core::cloning::ligation::tailing_recipe::Action::AddTail {
            base: TailingBase::A,
        };
        state.insert.treatment = Treatment::Phosphorylate;
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
            collagene_core::cloning::assembly::recipe::Preparation::WholeLinear
        ));
        assert_eq!(
            recipe.fragments[0].action,
            collagene_core::cloning::ligation::tailing_recipe::Action::AddTail {
                base: TailingBase::T
            }
        );
        assert_eq!(
            recipe.fragments[1].treatment,
            Some(EndTreatment::Phosphorylate)
        );
        assert_eq!(target.sequence.as_str(), "AACCTAACCA");
    }
    #[test]
    fn independently_declared_circle_and_linear_terminal_intent_feed_actual_tailing() {
        use collagene_core::{
            cloning::assembly::{
                intent as assembly_intent, junction_edit as assembly_junction_edit,
            },
            verification::evidence::Status,
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
            state.vector.action =
                collagene_core::cloning::ligation::tailing_recipe::Action::AddTail {
                    base: TailingBase::T,
                };
            state.insert.action =
                collagene_core::cloning::ligation::tailing_recipe::Action::AddTail {
                    base: TailingBase::A,
                };
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
                material_bindings: Vec::new(),
                name: "Operator declared target".into(),
                circular,
                fragments: roles
                    .iter()
                    .map(|role| collagene_core::cloning::assembly::intent::Contribution {
                        role: role.role.clone(),
                        input_id: role.input_id.clone(),
                        selector: collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                            start: 0,
                            end: 6,
                            reverse: false,
                        },
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
            let intent =
                assembly_junction_edit::build_intent(&inputs, &base, &edits, None).unwrap();
            assert_eq!(intent.target.sequence.as_str(), "AAACCCTGGGTTTA");
            let (bound, recipe) = tailing_request(&state, originals, &intent.target).unwrap();
            let preview = collagene_core::cloning::ligation::tailing_recipe::preview(
                &bound,
                &intent.target,
                &recipe,
            )
            .unwrap();
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
    fn kept_pre_tailed_vector_and_added_insert_use_one_recipe_without_double_tailing() {
        use collagene_core::{
            chemistry::Terminus,
            cloning::ligation::tailing_recipe::{self, Action, Profile, ProfileStatus},
            planning::plan::Operation,
            verification::evidence::Status,
        };
        let mut state = MethodsState::default();
        state.mode = Mode::Tailing;
        state.tailing_profile = Profile::Ta;
        state.vector.action = Action::KeepExisting;
        state.insert.action = Action::AddTail {
            base: TailingBase::A,
        };
        state.insert.treatment = Treatment::Phosphorylate;
        state.rows[0].preparation = OverlapPreparation::AnnealedOligos;
        state.rows[0].anneal_bottom_document = Some(2);
        state.rows[0].anneal_options.min_paired_bases = 4;
        state.rows[0].anneal_options.top_chemistry = Terminus::CLEAVED;
        state.rows[0].anneal_options.bottom_chemistry = Terminus::CLEAVED;
        state.rows[0].anneal_candidate = Some(0);
        state.rows[1].document = 1;
        state.rows[1].preparation = OverlapPreparation::ExactPcr;
        state.rows[1].forward = "GGAC".into();
        state.rows[1].reverse_primer = "GCAG".into();
        // Independently literal oligos establish 3′ T ends around paired AACC.
        // The insert starts as an 8nt blunt PCR duplex with disjoint 4nt
        // annealing domains; it receives the sole addition.
        let sources = vec![
            Molecule::dna("T-vector top", "AACCT", Topology::Linear).unwrap(),
            Molecule::dna("insert", "GGACCTGC", Topology::Linear).unwrap(),
            Molecule::dna("T-vector bottom", "GGTTT", Topology::Linear).unwrap(),
        ];
        let before = serde_json::to_vec(&sources).unwrap();
        let target =
            Molecule::dna("Independent circle", "AACCTGGACCTGCA", Topology::Circular).unwrap();
        let (inputs, recipe) = tailing_request(&state, sources.clone(), &target).unwrap();
        assert_eq!(inputs.len(), 3, "paired-oligo source remains present");
        assert_eq!(recipe.fragments[0].action, Action::KeepExisting);
        let revision = tailing_recipe::preview_revision(&inputs, &target, &recipe).unwrap();
        let preview = tailing_recipe::preview(&inputs, &target, &recipe).unwrap();
        preview
            .validate_for(&inputs, &target, &recipe, &revision)
            .unwrap();
        assert_eq!(
            preview.verification.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            preview.verification.findings
        );
        assert_eq!(
            preview.profile_status,
            ProfileStatus::Matched,
            "{:?}",
            preview.profile_findings
        );

        assert_eq!(preview.compiled.rows[0].tail_step, None);
        assert!(preview.compiled.rows[1].tail_step.is_some());
        assert_eq!(
            preview
                .compiled
                .plan
                .steps
                .iter()
                .filter(|step| matches!(step, Operation::Tail3Prime { .. }))
                .count(),
            1
        );
        assert_eq!(
            preview
                .verification
                .products
                .last()
                .unwrap()
                .sequence
                .as_str(),
            "AACCTGGACCTGCA"
        );
        assert_eq!(serde_json::to_vec(&sources).unwrap(), before);
        let mut changed = recipe.clone();
        changed.fragments[0].action = Action::AddTail {
            base: TailingBase::T,
        };
        assert!(
            preview
                .validate_for(&inputs, &target, &changed, &revision)
                .is_err()
        );
    }
    #[test]
    fn keeping_blunt_ends_can_ligate_without_satisfying_ta_profile() {
        use collagene_core::{
            cloning::ligation::tailing_recipe::{self, Action, Profile, ProfileStatus},
            verification::evidence::Status,
        };
        let mut state = MethodsState::default();
        state.mode = Mode::Tailing;
        state.tailing_profile = Profile::Ta;
        state.vector.action = Action::KeepExisting;
        state.insert.action = Action::KeepExisting;
        state.vector.treatment = Treatment::Phosphorylate;
        state.insert.treatment = Treatment::Phosphorylate;
        state.rows[0].preparation = OverlapPreparation::ExactPcr;
        state.rows[0].forward = "AAAC".into();
        state.rows[0].reverse_primer = "GCAG".into();
        state.rows[1].document = 1;
        state.rows[1].preparation = OverlapPreparation::ExactPcr;
        state.rows[1].forward = "GGGT".into();
        state.rows[1].reverse_primer = "TCGT".into();
        let inputs = vec![
            Molecule::dna("backbone", "AAACCTGC", Topology::Linear).unwrap(),
            Molecule::dna("cargo", "GGGTACGA", Topology::Linear).unwrap(),
        ];
        let target = Molecule::dna(
            "Independent blunt circle",
            "AAACCTGCGGGTACGA",
            Topology::Circular,
        )
        .unwrap();
        let (inputs, recipe) = tailing_request(&state, inputs, &target).unwrap();
        let preview = tailing_recipe::preview(&inputs, &target, &recipe).unwrap();
        assert!(
            preview
                .compiled
                .rows
                .iter()
                .all(|row| row.tail_step.is_none())
        );
        assert_eq!(
            preview.verification.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            preview.verification.findings
        );
        assert_eq!(preview.profile_status, ProfileStatus::Incompatible);

        assert_eq!(
            preview
                .verification
                .products
                .last()
                .unwrap()
                .sequence
                .as_str(),
            "AAACCTGCGGGTACGA"
        );
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
    use crate::cloning::preparation::cloning_methods_model::Mode;
    use crate::cloning::preparation::cloning_methods_model::Treatment;
    use collagene_core::{
        cloning::{assembly::intent as assembly_intent, restriction::catalog as catalog_digest},
        verification::evidence::Status,
    };
    use collagene_model::Topology;
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
            material_bindings: Vec::new(),
            name: "Independent backbone and top contribution".into(),
            circular: true,
            fragments: vec![
                collagene_core::cloning::assembly::intent::Contribution {
                    role: roles[0].role.clone(),
                    input_id: ids[0].clone(),
                    selector: collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                        start: 2,
                        end: 14,
                        reverse: false,
                    },
                },
                collagene_core::cloning::assembly::intent::Contribution {
                    role: roles[1].role.clone(),
                    input_id: ids[1].clone(),
                    selector: collagene_core::cloning::assembly::intent::Selector::ReferenceRange {
                        start: 0,
                        end: 10,
                        reverse: false,
                    },
                },
            ],
        };
        let intent = assembly_intent::build(&inputs, &declaration, None).unwrap();
        assert_eq!(intent.target.sequence.as_str(), "AATTCTTACGGGAATTGCTACG");
        let (bound_ids, fragments) = prepared_rows(&state, 3).unwrap();
        assert_eq!(bound_ids, ids);
        assert_eq!(fragments.len(), 2);
        assert!(
            matches!(&fragments[1].preparation,collagene_core::cloning::assembly::recipe::Preparation::AnnealedOligos{bottom_input_id,..} if bottom_input_id==&ids[2])
        );
        let recipe = collagene_core::cloning::ligation::recipe::Recipe {
            name: "Actual three-source ligation".into(),
            fragments: fragments
                .into_iter()
                .zip(&state.rows)
                .map(
                    |(source, row)| collagene_core::cloning::ligation::recipe::Fragment {
                        source,
                        treatment: row.treatment.operation(),
                    },
                )
                .collect(),
            circular: true,
            allow_nicked: false,
        };
        let preview = collagene_core::cloning::ligation::recipe::preview(
            &inputs,
            &intent.target,
            &recipe,
            collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
        )
        .unwrap();
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
        state.source_roles = Some(vec![crate::assembly::draft::assembly_role_ids::Role {
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
        let collagene_core::cloning::assembly::recipe::Preparation::MismatchPcr {
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
            collagene_core::cloning::pcr::mutagenesis::PcrOptions::default().max_pair_checks
        );
        assert!(state == before);
        state.rows[0].mismatch_max_mismatches = "99".into();
        assert!(prepared_rows(&state, 1).is_err());
    }
}

#[cfg(test)]
mod overlap_extension_tests {
    use super::*;
    use crate::cloning::preparation::cloning_methods_model::ExtensionPrimerInput;
    use crate::cloning::preparation::cloning_methods_model::Mode;
    use crate::cloning::preparation::cloning_methods_model::OverlapExtensionPcrInput;
    use collagene_core::{
        cloning::assembly::recipe as assembly_recipe,
        documents::workspace::Workspace,
        material::{
            lineage::{EventBinding, InputEventBinding},
            outputs as cloning_outputs, run as run_material,
        },
        planning::plan::{DesignPlan, Operation, Source},
        verification::evidence::Status,
    };
    use collagene_model::Topology;
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
            crate::assembly::draft::assembly_role_ids::Role {
                role: "edited-insert".into(),
                input_id: "original-stock-event".into(),
            },
            crate::assembly::draft::assembly_role_ids::Role {
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
                template_strand: collagene_core::chemistry::Strand::Top,
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
        let capsule =
            collagene_core::material::capsule::capsule(&commit.workspace, &selection).unwrap();
        let saved = commit.workspace.revisions[&selection.event.revision].clone();
        let mut s = state();
        s.rows[0].run_origin = Some(
            crate::history::run_origin_model::PinnedOrigin::from_material(
                collagene_core::material::run::MaterialOrigin::Run(capsule.clone()),
            )
            .unwrap(),
        );
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
        let catalog = collagene_core::cloning::restriction::catalog::builtin_catalog().unwrap();
        let id = catalog
            .edition()
            .enzymes
            .iter()
            .find(|entry| entry.name == "EcoRI")
            .unwrap()
            .id
            .clone();
        s.rows[0].catalog_selection =
            Some(collagene_core::cloning::restriction::catalog::select(&catalog, &[id]).unwrap());
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
        let catalog = collagene_core::cloning::restriction::catalog::builtin_catalog().unwrap();
        let ids = catalog
            .edition()
            .enzymes
            .iter()
            .filter(|entry| entry.name == "EcoRI" || entry.name == "XhoI")
            .map(|entry| entry.id.clone())
            .collect::<Vec<_>>();
        row.catalog_selection =
            Some(collagene_core::cloning::restriction::catalog::select(&catalog, &ids).unwrap());
        let before = s.clone();
        let (_, rows) = prepared_rows(&s, 2).unwrap();
        let assembly_recipe::Preparation::OverlapExtensionPcrCatalogDigest {
            conditions: _,
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

#[cfg(test)]
mod source_generation_tests {
    use super::*;
    use crate::cloning::preparation::cloning_methods_model::SourceGenerationInput;
    use collagene_core::{
        chemistry::tracked_duplex::ModificationSummary,
        planning::plan::{DesignPlan, Operation, PlanInput},
        verification,
    };
    use collagene_model::{Molecule, Topology};
    #[test]
    fn original_generation_precedes_pcr_and_new_strands_do_not_copy_modifications() {
        let source = Molecule::dna(
            "Dam declared original",
            "CAGTGATCGTACGTAC",
            Topology::Linear,
        )
        .unwrap();
        let target = Molecule::dna(
            "Independent intended PCR",
            "CAGTGATCGTACGTAC",
            Topology::Linear,
        )
        .unwrap();
        let mut state = MethodsState::default();
        state.rows.truncate(1);
        let row = &mut state.rows[0];
        row.preparation = OverlapPreparation::ExactPcr;
        row.forward = "CAGTGATC".into();
        row.reverse_primer = "GTACGTAC".into();
        row.source_generation = Some(SourceGenerationInput {
            dam: "fully".into(),
            dcm: "unmodified".into(),
            cpg: "unmodified".into(),
        });
        let (ids, fragments) = prepared_rows(&state, 1).unwrap();
        assert!(fragments[0].generation.is_some());
        let inputs = vec![PlanInput {
            id: ids[0].clone(),
            revision: source.revision(),
            molecule: source,
        }];
        let compiled = collagene_core::cloning::assembly::recipe::compile_preparations(
            &inputs,
            &target,
            "Original declaration before PCR",
            &fragments,
            collagene_core::cloning::assembly::recipe::CleavagePolicy::PositionalChecked,
        )
        .unwrap();
        assert!(matches!(
            &compiled.steps[0],
            Operation::DeclareSourceModifications { .. }
        ));
        assert!(matches!(&compiled.steps[1], Operation::Pcr { .. }));
        let plan = DesignPlan {
            name: "Actual generated strands".into(),
            target_revision: target.revision(),
            inputs,
            steps: compiled.steps,
        };
        let report = verification::engine::verify(&plan, &target);
        assert_eq!(report.products.len(), 2, "{:?}", report.findings);
        let before: ModificationSummary =
            serde_json::from_str(&report.products[0].metadata["tracked_modifications"]).unwrap();
        let after: ModificationSummary =
            serde_json::from_str(&report.products[1].metadata["tracked_modifications"]).unwrap();
        assert!(before.known_modified > 0);
        assert_eq!(after.known_modified, 0);
        assert_eq!(after.unknown_eligible, 0);
        assert_eq!(report.products[1].sequence.as_str(), "CAGTGATCGTACGTAC");
    }
}

pub(crate) fn overlap_request(
    state: &MethodsState,
    inputs: Vec<Molecule>,
    target: &Molecule,
) -> Result<
    (
        Vec<PlanInput>,
        collagene_core::cloning::assembly::recipe::Recipe,
    ),
    String,
> {
    use collagene_core::cloning::assembly::recipe::Recipe;
    state.validate_processing_selection()?;
    let (overlaps, options) = state.profile()?;
    let (ids, fragments) = prepared_rows(state, inputs.len())?;
    let recipe = Recipe {
        name: format!("Profiled overlap assembly to {}", target.name),
        fragments,
        overlaps,
        circular: state.circular,
        options,
    };
    let borrowed = ids
        .iter()
        .zip(&inputs)
        .map(|(id, m)| (id.as_str(), m))
        .collect::<Vec<_>>();
    if let Some(form) = state.annealing.as_ref().filter(|form| form.enabled) {
        let policy = form.policy(
            recipe.overlaps.clone(),
            recipe.circular,
            recipe.fragments.len(),
        )?;
        collagene_core::cloning::assembly::recipe::annealing_source_preflight(
            &borrowed, target, &recipe, &policy,
        )
        .map_err(|e| e.to_string())?;
    } else if let Some(form) = state.homology.as_ref().filter(|form| form.enabled) {
        let policy = form.policy(
            recipe.overlaps.clone(),
            recipe.circular,
            recipe.fragments.len(),
        )?;
        collagene_core::cloning::assembly::recipe::homology_source_preflight(
            &borrowed, target, &recipe, &policy,
        )
        .map_err(|e| e.to_string())?;
    } else {
        collagene_core::cloning::assembly::recipe::source_preflight(&borrowed, target, &recipe)
            .map_err(|e| e.to_string())?;
    }
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

pub(crate) fn ligation_request(
    state: &MethodsState,
    inputs: Vec<Molecule>,
    target: &Molecule,
) -> Result<
    (
        Vec<PlanInput>,
        collagene_core::cloning::ligation::recipe::Recipe,
    ),
    String,
> {
    use collagene_core::cloning::ligation::recipe::{Fragment, Recipe};
    let (ids, fragments) = prepared_rows(state, inputs.len())?;
    if state.mode == Mode::GoldenGate {
        let geometry = state.golden_gate.enzyme.profile();
        for fragment in &fragments {
            let selection=fragment.preparation.catalog_selection()
                .ok_or("Golden Gate rows require a pinned Type IIS catalog digest, including an explicit exact or intentional-mismatch PCR followed by digestion.")?;
            if selection.enzymes.len() != 1
                || !selection.enzymes.iter().all(|entry| {
                    entry.enzyme.recognition == geometry.recognition
                        && entry.enzyme.top_cut == geometry.top_cut
                        && entry.enzyme.bottom_cut == geometry.bottom_cut
                })
            {
                return Err("A row's pinned enzyme geometry differs from the chosen Type IIS enzyme. Choose the matching existing digest fragment or apply fresh PCR-site proposals.".into());
            }
        }
    }
    let recipe = Recipe {
        name: format!(
            "{} to {}",
            if state.mode == Mode::GoldenGate {
                "Ordered Type IIS restriction/ligation"
            } else {
                "Ordered ligation"
            },
            target.name
        ),
        fragments: fragments
            .into_iter()
            .zip(&state.rows)
            .map(|(source, row)| Fragment {
                source,
                treatment: row.treatment.operation(),
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
    collagene_core::cloning::ligation::recipe::source_preflight(&borrowed, target, &recipe)
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
