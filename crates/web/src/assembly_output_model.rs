//! Optional output capabilities from the existing native order extractor.
//! Collection availability is unrelated to successful molecular reconstruction.
use collagene_core::plan_oligos;
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct PrimerCapability {
    pub distinct_oligos: usize,
    pub uses: usize,
    pub complete: bool,
    pub unavailable: Vec<String>,
    pub warnings: Vec<String>,
}
impl PrimerCapability {
    pub(crate) fn from_report(
        report: &plan_oligos::Report,
        plan_revision: &str,
        target_revision: &str,
    ) -> Result<Self, String> {
        report.validate().map_err(|e| e.to_string())?;
        if report.plan_revision != plan_revision || report.target_revision != target_revision {
            return Err("Oligo capabilities belong to a different plan or intended target.".into());
        }
        Ok(Self {
            distinct_oligos: report.oligos.records.len(),
            uses: report.uses.len(),
            complete: report.complete,
            unavailable: report.unavailable.clone(),
            warnings: report.warnings.clone(),
        })
    }
    pub(crate) fn can_save(&self) -> bool {
        self.complete && self.distinct_oligos > 0
    }
    pub(crate) fn selected(&self, requested: bool, context_current: bool) -> bool {
        requested && context_current && self.can_save()
    }
    pub(crate) fn summary(&self) -> String {
        if !self.complete {
            format!(
                "Primer collection unavailable: {} declared oligos were found, but this plan's order extraction is incomplete. The construct can still be created without a primer collection.",
                self.distinct_oligos
            )
        } else if self.distinct_oligos == 0 {
            "No new orderable oligos in this recipe. Digest-only preparations and already-prepared material do not imply a new primer order. Create the construct without a primer collection.".into()
        } else {
            format!(
                "{} distinct declared oligos · {} recorded uses. Saving a collection does not establish synthesis or experimental testing.",
                self.distinct_oligos, self.uses
            )
        }
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::{
        Molecule, Topology,
        analysis::Primer,
        cloning_outputs,
        plan_oligos::Options,
        run_lineage::{EventBinding, InputEventBinding},
        run_material,
        verification::{self, DesignPlan, Operation, PlanInput, Source, TopoSource},
        workspace::Workspace,
    };
    fn fixture() -> (DesignPlan, Molecule) {
        let source = Molecule::dna("Original16bp", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
        let target =
            Molecule::dna("Independent16bp", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
        let plan = DesignPlan {
            name: "Declared actual primers".into(),
            target_revision: target.revision(),
            inputs: vec![PlanInput {
                id: "original".into(),
                revision: source.revision(),
                molecule: source,
            }],
            steps: vec![Operation::Pcr {
                template_strand: collagene_core::chemistry::Strand::Top,
                template: Source::Input("original".into()),
                forward: Primer::new("F", "ACGACCCC", "").unwrap(),
                reverse: Primer::new("R", "GTTACCCC", "").unwrap(),
            }],
        };
        (plan, target)
    }
    fn capability(plan: &DesignPlan, target: &Molecule) -> PrimerCapability {
        let report = plan_oligos::build(plan, target, &Options::default()).unwrap();
        PrimerCapability::from_report(
            &report,
            &verification::plan_revision(plan),
            &target.revision(),
        )
        .unwrap()
    }
    #[test]
    fn real_empty_and_nonempty_native_orders_control_only_the_optional_collection() {
        let (mut plan, target) = fixture();
        let available = capability(&plan, &target);
        assert_eq!((available.distinct_oligos, available.uses), (2, 2));
        assert!(available.can_save() && available.selected(true, true));
        assert!(!available.selected(true, false) && !available.selected(false, true));
        plan.steps = vec![Operation::Digest {
            template: Source::Input("original".into()),
            enzymes: vec!["SmaI".into()],
            fragment_index: 0,
        }];
        let empty = capability(&plan, &target);
        assert_eq!(empty.distinct_oligos, 0);
        assert!(!empty.can_save() && !empty.selected(true, true));
        assert!(empty.summary().contains("without a primer collection"));
    }
    #[test]
    fn real_embedded_topo_primers_enable_collection_without_verifying_the_graph() {
        let (mut plan, target) = fixture();
        let Operation::Pcr {
            forward, reverse, ..
        } = plan.steps[0].clone()
        else {
            unreachable!()
        };
        plan.steps.push(Operation::Topo {
            sources: vec![TopoSource {
                role: "template".into(),
                source: Source::Input("original".into()),
            }],
            request: collagene_core::topo::Request {
                recipe: collagene_core::topo::GenerationRecipe {
                    schema_version: 1,
                    materials: vec![collagene_core::topo::MaterialRecipe::ExactPcr {
                        material_id: "embedded-PCR".into(),
                        template_role: "template".into(),
                        forward,
                        reverse,
                        tail: None,
                    }],
                },
                transfers: vec![],
                engine_limits: Default::default(),
                projection_options: Default::default(),
            },
            selection: collagene_core::topo::Selection {
                batch_revision: "0".repeat(64),
                component: 0,
                reference_seed: collagene_core::topo::physical::Call {
                    strand: 0,
                    index: 0,
                },
            },
        });
        let cap = capability(&plan, &target);
        assert_eq!(cap.distinct_oligos, 2);
        assert_eq!(cap.uses, 4);
        assert!(cap.complete && cap.can_save() && cap.selected(true, true));
        assert!(cap.unavailable.is_empty());
        assert!(!cap.selected(true, false));
        // Literal primer declarations do not validate this deliberately unbound graph.
        assert_ne!(
            collagene_core::verification::verify(&plan, &target).status,
            collagene_core::verification::Status::ConsistentUnderAssumptions
        );
    }
    fn binding(w: &Workspace, i: usize) -> EventBinding {
        let d = &w.documents[i];
        EventBinding {
            document_id: d.id.clone(),
            event_id: d.head_event.clone().unwrap(),
            revision: d.head.clone(),
        }
    }
    #[test]
    fn historical_origin_primers_are_excluded_and_only_downstream_new_primers_enable_collection() {
        let (plan, target) = fixture();
        let w = Workspace::from_documents(
            "Original and intent",
            vec![plan.inputs[0].molecule.clone(), target.clone()],
            1,
        )
        .unwrap();
        let commit = cloning_outputs::finalize(
            &w,
            &plan,
            binding(&w, 1),
            vec![InputEventBinding {
                input_id: "original".into(),
                event: binding(&w, 0),
            }],
            &cloning_outputs::Options {
                name: "Named saved material".into(),
                show_prepared: false,
                save_primers: false,
                allow_incomplete: false,
            },
        )
        .unwrap();
        let event = commit.receipt.construct;
        let origin = run_material::capsule(
            &commit.workspace,
            &run_material::Selection {
                run_id: commit.receipt.run.run_id,
                step_index: 0,
                event: event.clone(),
            },
        )
        .unwrap();
        let molecule = commit.workspace.revisions[&event.revision].clone();
        let mut next = DesignPlan {
            name: "Reuse already-prepared material".into(),
            target_revision: target.revision(),
            inputs: vec![PlanInput {
                id: "saved".into(),
                revision: molecule.revision(),
                molecule,
            }],
            steps: vec![Operation::ReconstructRunProduct {
                template: Source::Input("saved".into()),
                origin,
            }],
        };
        let cap = capability(&next, &target);
        assert_eq!(cap.distinct_oligos, 0);
        assert!(!cap.can_save());
        assert!(
            cap.warnings
                .iter()
                .any(|s| s.contains("Historical origin PCR primers"))
        );
        let Operation::Pcr {
            forward, reverse, ..
        } = plan.steps[0].clone()
        else {
            unreachable!()
        };
        next.steps.push(Operation::Pcr {
            template_strand: collagene_core::chemistry::Strand::Top,
            template: Source::Step(0),
            forward,
            reverse,
        });
        let cap = capability(&next, &target);
        assert_eq!((cap.distinct_oligos, cap.uses), (2, 2));
        assert!(cap.can_save());
        let report = plan_oligos::build(&next, &target, &Options::default()).unwrap();
        assert!(
            PrimerCapability::from_report(
                &report,
                &verification::plan_revision(&plan),
                &target.revision()
            )
            .is_err()
        );
    }
}
