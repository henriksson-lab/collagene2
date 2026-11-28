//! Bounded freshly evaluated documentary evidence across material ancestry.
//! These transport validators establish historical integrity, never current science.
use crate::{
    cloning::restriction::conditions::source::SourceReport,
    verification::evidence::{PhysicalProcessingEvidence, VerificationReport},
};
use collagene_model::Error;
use std::collections::BTreeMap;

pub const MAX_CONTEXTS: usize = 32;
pub const MAX_BYTES: usize = 4 * 1024 * 1024;

pub fn validate(reports: &[SourceReport]) -> Result<(), Error> {
    if reports.len() > MAX_CONTEXTS {
        return Err(Error::TooLarge);
    }
    crate::encoding::bounded_json::preflight(&reports, MAX_BYTES)?;
    let mut contexts = std::collections::BTreeSet::new();
    for report in reports {
        report.validate_historical()?;
        if !contexts.insert(&report.context_revision) {
            return Err(super::bad("duplicate documentary condition context"));
        }
    }
    Ok(())
}
/// Resource checks precede cloning/combining full records. Identical contexts
/// shared by multiple ancestry paths are retained once; inconsistent duplicates fail.
pub(crate) fn collect<'a>(
    reports: impl IntoIterator<Item = &'a SourceReport>,
) -> Result<Vec<SourceReport>, Error> {
    let mut distinct = BTreeMap::<&str, &SourceReport>::new();
    for report in reports {
        if let Some(previous) = distinct.get(report.context_revision.as_str()) {
            if *previous != report {
                return Err(super::bad("conflicting documentary condition context"));
            }
        } else {
            if distinct.len() >= MAX_CONTEXTS {
                return Err(Error::TooLarge);
            }
            distinct.insert(&report.context_revision, report);
        }
    }
    let borrowed = distinct.into_values().collect::<Vec<_>>();
    crate::encoding::bounded_json::preflight(&borrowed, MAX_BYTES)?;
    for report in &borrowed {
        report.validate_historical()?;
    }
    Ok(borrowed.into_iter().cloned().collect())
}
/// Collect documentary evidence only along the selected operation's typed source
/// graph. The capsule binds the current report to its exact candidate and target;
/// neither this inspection nor producer reports create material/verification authority.
pub fn from_origin_report(
    origin: &crate::material::capsule::Capsule,
    report: &VerificationReport,
) -> Result<Vec<SourceReport>, Error> {
    fn walk(
        origin: &crate::material::capsule::Capsule,
        report: &VerificationReport,
        depth: usize,
        nodes: &mut usize,
        bytes: &mut usize,
    ) -> Result<Vec<SourceReport>, Error> {
        use crate::planning::plan::Source;
        *nodes = nodes.checked_add(1).ok_or(Error::TooLarge)?;
        *bytes = bytes
            .checked_add(origin.candidate.len())
            .filter(|n| *n <= super::MAX_CAPSULE_BYTES)
            .ok_or(Error::TooLarge)?;
        if depth > 32 || *nodes > 128 {
            return Err(Error::TooLarge);
        }
        // Current envelope/closed-origin admission precedes decoding candidate
        // collections. This does not execute the plan or trust its claimed result.
        origin.validate_documentary()?;
        report.validate_documentary()?;
        let candidate = crate::documents::transfer::decode_candidate_shallow(&origin.candidate)?;
        let selected = origin.selection.step_index;
        if report.engine_version != crate::verification::ENGINE_VERSION
            || report.plan_revision != candidate.plan_revision
            || candidate.plan_revision != candidate.plan.revision()
            || report.target_revision != candidate.plan.target_revision
            || report.target_revision != origin.target_event.revision
            || report.target_revision != origin.target.revision()
            || report.status == crate::verification::evidence::Status::Invalid
            || report.steps.len() != origin.products.len()
            || selected >= candidate.plan.steps.len()
            || !report
                .steps
                .iter()
                .zip(&origin.products)
                .all(|(step, binding)| {
                    step.state == crate::verification::evidence::StepState::Reconstructed
                        && step
                            .product
                            .as_ref()
                            .is_some_and(|product| product.revision() == binding.revision)
                })
        {
            return Err(super::bad(
                "documentary conditions differ from the exact selected origin plan/target/product",
            ));
        }
        let mut reachable = std::collections::BTreeSet::from([selected]);
        let mut pending = vec![selected];
        while let Some(index) = pending.pop() {
            let operation = candidate.plan.steps.get(index).ok_or_else(|| {
                super::bad("selected condition ancestry references an absent step")
            })?;
            for source in operation.sources() {
                match source {
                    Source::Step(parent) => {
                        if *parent >= index {
                            return Err(super::bad(
                                "condition ancestry requires an earlier source step",
                            ));
                        }
                        if reachable.insert(*parent) {
                            pending.push(*parent);
                        }
                    }
                    Source::Input(id) => {
                        if !candidate.plan.inputs.iter().any(|input| input.id == *id) {
                            return Err(super::bad(
                                "condition ancestry references an absent input",
                            ));
                        }
                    }
                }
            }
        }
        let mut found = collect(
            report
                .step_reports
                .iter()
                .filter(|step| reachable.contains(&step.step_index))
                .filter_map(|step| step.reaction_conditions.as_ref()),
        )?;
        for step in report
            .physical_processing_reports
            .iter()
            .filter(|step| reachable.contains(&step.step_index))
        {
            match &step.evidence {
                PhysicalProcessingEvidence::DerivedOrigin(evidence) => {
                    found = collect(found.iter().chain(&evidence.reaction_conditions))?;
                }
                PhysicalProcessingEvidence::RunProduct(evidence) => {
                    *bytes = bytes
                        .checked_add(evidence.current_report.len())
                        .filter(|n| *n <= super::MAX_CAPSULE_BYTES)
                        .ok_or(Error::TooLarge)?;
                    let nested = evidence.decode_report()?;
                    let inherited = walk(&evidence.origin, &nested, depth + 1, nodes, bytes)?;
                    found = collect(found.iter().chain(&inherited))?;
                }
                _ => {}
            }
        }
        Ok(found)
    }
    crate::encoding::bounded_json::preflight(report, super::MAX_CAPSULE_BYTES)?;
    let mut nodes = 0;
    let mut bytes = 0;
    walk(origin, report, 0, &mut nodes, &mut bytes)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        chemistry::{Duplex, tracked_duplex::TrackedDuplex},
        cloning::restriction::{
            catalog,
            conditions::{ChosenFormulation, ReactionConditions},
        },
    };
    use collagene_model::{Molecule, Topology};
    fn report(minutes: u32) -> SourceReport {
        let source = Molecule::dna("original", "AAGAATTCTT", Topology::Linear).unwrap();
        let material =
            TrackedDuplex::from_unknown_duplex(&Duplex::from_molecule(&source).unwrap()).unwrap();
        let catalog = catalog::builtin_catalog().unwrap();
        let selection = catalog::select(&catalog, &["ecori".into()]).unwrap();
        let options = ReactionConditions {
            schema_version: 1,
            formulations: vec![ChosenFormulation {
                enzyme_id: "ecori".into(),
                formulation_id: Some("neb-r3101".into()),
            }],
            buffer: "rCutSmart".into(),
            temperature_milli_c: 37000,
            incubation_minutes: minutes,
            minimum_activity_percent: 75,
            activity_context: Some("NEB unit-assay relative activity".into()),
            temperature_context: Some("vendor recommended incubation".into()),
        };
        SourceReport::assess(&catalog, &source, &material, &selection, &options).unwrap()
    }
    #[test]
    fn repeated_exact_context_deduplicates_but_inconsistent_duplicate_is_rejected() {
        let original = report(60);
        assert_eq!(
            collect([&original, &original]).unwrap(),
            vec![original.clone()]
        );
        let mut forged = original.clone();
        forged.material_revision = "0".repeat(64);
        assert!(collect([&original, &forged]).is_err());
    }
    #[test]
    fn distinct_context_and_serialized_byte_limits_refuse_without_truncation() {
        let reports = (1..=33).map(report).collect::<Vec<_>>();
        assert_eq!(collect(reports[..32].iter()).unwrap().len(), 32);
        assert!(matches!(collect(reports.iter()), Err(Error::TooLarge)));
        let mut giant = report(60);
        giant.source_revision = "a".repeat(MAX_BYTES);
        assert!(matches!(collect([&giant]), Err(Error::TooLarge)));
    }
    /// Documentary capsule assertions use shaped test events; actual product
    /// sequences and evidence come from the ordinary fresh deterministic verifier.
    fn branching_origin(shared: bool) -> (crate::material::capsule::Capsule, VerificationReport) {
        use crate::{
            chemistry::EndRepairMethod,
            documents::transfer::{Candidate, save_candidate},
            material::{
                lineage::{EventBinding, InputEventBinding, ProductBinding},
                run::Selection,
            },
            planning::plan::{DesignPlan, Operation, PlanInput, Source},
        };
        let source = Molecule::dna("source", "AAGAATTCAAACAATTC", Topology::Linear).unwrap();
        let target = Molecule::dna(
            "independent target",
            if shared {
                "AATTCAAACAATTCAAACAATTC"
            } else {
                "AATTCAAACAATTC"
            },
            Topology::Linear,
        )
        .unwrap();
        let catalog = catalog::builtin_catalog().unwrap();
        let selection = catalog::select(&catalog, &["ecori".into()]).unwrap();
        let options = |minutes| ReactionConditions {
            incubation_minutes: minutes,
            ..report(60).report.options
        };
        let mut steps = vec![
            Operation::CatalogDigest {
                template: Source::Input("s".into()),
                selection: selection.clone(),
                conditions: Some(options(99)),
                fragment_index: 1,
            },
            Operation::CatalogDigest {
                template: Source::Input("s".into()),
                selection,
                conditions: Some(options(60)),
                fragment_index: 1,
            },
            Operation::RepairEnds {
                template: Source::Step(1),
                method: EndRepairMethod::FillFivePrimeOverhangs,
            },
        ];
        if shared {
            steps.push(Operation::RepairEnds {
                template: Source::Step(1),
                method: EndRepairMethod::FillFivePrimeOverhangs,
            });
            steps.push(Operation::Overlap {
                fragments: vec![Source::Step(2), Source::Step(3)],
                overlaps: vec![5],
                circular: false,
            });
        }
        let plan = DesignPlan {
            name: "Disconnected digest plus selected dependency path".into(),
            target_revision: target.revision(),
            inputs: vec![PlanInput {
                id: "s".into(),
                revision: source.revision(),
                molecule: source.clone(),
            }],
            steps,
        };
        let report = crate::verification::engine::verify(&plan, &target);
        assert_eq!(
            report.steps.len(),
            plan.steps.len(),
            "{:?}",
            report.findings
        );
        assert_ne!(
            report.status,
            crate::verification::evidence::Status::Invalid,
            "{:?}",
            report.findings
        );
        let selected = plan.steps.len() - 1;
        let products = report
            .steps
            .iter()
            .map(|step| {
                let index = step.step_index;
                let product = step.product.as_ref().unwrap();
                ProductBinding {
                    step_index: index,
                    document_id: format!("product-{index}"),
                    event_id: format!("event-{:064x}", index + 1),
                    revision: product.revision(),
                }
            })
            .collect::<Vec<_>>();
        let product = &products[selected];
        let origin = crate::material::capsule::Capsule {
            schema_version: 1,
            selection: Selection {
                run_id: "a".repeat(64),
                step_index: selected,
                event: EventBinding {
                    document_id: product.document_id.clone(),
                    event_id: product.event_id.clone(),
                    revision: product.revision.clone(),
                },
            },
            candidate: save_candidate(&Candidate::new(plan)).unwrap(),
            target_event: EventBinding {
                document_id: "target".into(),
                event_id: format!("event-{:064x}", 100),
                revision: target.revision(),
            },
            target,
            inputs: vec![InputEventBinding {
                input_id: "s".into(),
                event: EventBinding {
                    document_id: "source".into(),
                    event_id: format!("event-{:064x}", 101),
                    revision: source.revision(),
                },
            }],
            products,
        };
        (origin, report)
    }
    #[test]
    fn selected_dependency_path_excludes_disconnected_earlier_digest() {
        let (origin, report) = branching_origin(false);
        assert_eq!(
            report
                .step_reports
                .iter()
                .filter(|s| s.reaction_conditions.is_some())
                .count(),
            2
        );
        let selected = from_origin_report(&origin, &report).unwrap();
        assert_eq!(selected.len(), 1);
        assert_eq!(selected[0].report.options.incubation_minutes, 60);
        assert!(
            !selected
                .iter()
                .any(|e| e.report.options.incubation_minutes == 99)
        );
        let mut wrong = report.clone();
        wrong.plan_revision = "0".repeat(64);
        assert!(from_origin_report(&origin, &wrong).is_err());
        wrong = report.clone();
        wrong.target_revision = "0".repeat(64);
        assert!(from_origin_report(&origin, &wrong).is_err());
        wrong = report.clone();
        wrong.steps[origin.selection.step_index]
            .product
            .as_mut()
            .unwrap()
            .name = "Tampered verifier output".into();
        assert!(from_origin_report(&origin, &wrong).is_err());
        let mut another = origin;
        another.selection.step_index = another.products.len();
        assert!(from_origin_report(&another, &report).is_err());
    }
    #[test]
    fn renamed_outer_event_does_not_replace_recorded_verifier_product_identity() {
        let (mut origin, report) = branching_origin(false);
        let mut renamed = report.product(origin.selection.step_index).unwrap().clone();
        renamed.name = "Saved actual source".into();
        assert_ne!(
            renamed.revision(),
            origin.products[origin.selection.step_index].revision
        );
        origin.selection.event.revision = renamed.revision();
        let retained = from_origin_report(&origin, &report).unwrap();
        assert_eq!(retained.len(), 1);
        assert_eq!(retained[0].report.options.incubation_minutes, 60);
    }
    #[test]
    fn shared_selected_ancestor_is_retained_once_and_unrelated_branch_stays_excluded() {
        let (origin, report) = branching_origin(true);
        let selected = from_origin_report(&origin, &report).unwrap();
        assert_eq!(selected.len(), 1);
        assert_eq!(selected[0].report.options.incubation_minutes, 60);
    }
}
