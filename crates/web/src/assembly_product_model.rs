//! Bounded readonly views of a host-supplied fresh reconstruction snapshot.
//! Prefix documents never grant final-product, installation or material authority.
use collagene_core::{
    chemistry::DuplexSummary,
    verification::{self, DesignPlan, Operation, Source, Status, VerificationReport},
};
use collagene_model::Molecule;
fn clip(value: &str, limit: usize) -> String {
    let mut chars = value.chars();
    let text = chars.by_ref().take(limit).collect::<String>();
    if chars.next().is_some() {
        format!("{text}…")
    } else {
        text
    }
}
fn source_label(source: &Source) -> String {
    match source {
        Source::Input(id) => format!("original input {}", clip(id, 80)),
        Source::Step(index) => format!("step {}", index + 1),
    }
}
pub(crate) fn operation_label(operation: &Operation) -> String {
    match operation {
        Operation::OverlapExtension {
            left,
            right,
            options,
        } => format!(
            "Reciprocal overlap extension of {} / {} · {} bp terminal pairing",
            source_label(&left.source),
            source_label(&right.source),
            options.overlap
        ),
        Operation::Pcr {
            template_strand: _,
            template,
            forward,
            reverse,
        } => format!(
            "Exact PCR of {} · {} / {}",
            source_label(template),
            clip(&forward.name, 60),
            clip(&reverse.name, 60)
        ),
        Operation::PcrWithMismatches {
            template,
            forward,
            reverse,
            options,
        } => format!(
            "Intentional-mismatch PCR of {} · {} / {} · at most {} mismatches, {} exact 3′ bases",
            source_label(template),
            clip(&forward.name, 60),
            clip(&reverse.name, 60),
            options.binding.max_mismatches,
            options.binding.min_three_prime_exact
        ),
        Operation::Digest {
            template,
            enzymes,
            fragment_index,
        } => format!(
            "Digest {} · {} · fragment {} (zero based)",
            source_label(template),
            enzymes
                .iter()
                .take(8)
                .map(|name| clip(name, 32))
                .collect::<Vec<_>>()
                .join(", "),
            fragment_index
        ),
        Operation::CatalogDigest {
            template,
            fragment_index,
            ..
        } => format!(
            "Catalog-pinned digest of {} · fragment {} (zero based)",
            source_label(template),
            fragment_index
        ),
        Operation::RepairEnds { template, method } => format!(
            "Explicit end repair of {} · {:?}",
            source_label(template),
            method
        ),
        Operation::TreatEnds {
            template,
            treatment,
        } => format!(
            "Explicit terminal treatment of {} · {:?}",
            source_label(template),
            treatment
        ),
        Operation::DeclareSourceTermini {
            template,
            left,
            right,
        } => format!(
            "Original duplex terminal assumption for {} · left {:?} · right {:?}",
            source_label(template),
            left,
            right
        ),
        Operation::AnnealOligos {
            top,
            bottom,
            candidate_index,
            ..
        } => format!(
            "Anneal ordered oligos {} / {} · candidate {} (zero based)",
            source_label(top),
            source_label(bottom),
            candidate_index
        ),
        Operation::Tail3Prime { template, base } => format!(
            "Add explicit {} 3′ tails to {}",
            base.as_str(),
            source_label(template)
        ),
        Operation::ReconstructMaterial { template, origin } => format!(
            "Reconstruct reusable PCR/digest material from {} · exact event {}",
            source_label(template),
            origin.product.event_id
        ),

        Operation::ReconstructRunProduct { template, origin } => format!(
            "Reconstruct saved material {} · origin run {} · product {}",
            source_label(template),
            clip(&origin.selection.run_id, 16),
            origin.selection.step_index + 1
        ),
        Operation::HomologyAssembly {
            fragments, policy, ..
        } => format!(
            "Reconstruct homology-assembly strands · {} fragments · {} · joining homology {:?} bp",
            fragments.len(),
            if policy.topology == collagene_model::Topology::Circular {
                "circular"
            } else {
                "linear"
            },
            policy.overlaps
        ),
        Operation::ResectThreePrime { cuts, .. } => format!(
            "Explicit 3′ resection · left bottom suffix {} bp · right top suffix {} bp",
            cuts.left, cuts.right
        ),
        Operation::HomologyAnnealing {
            fragments, policy, ..
        } => format!(
            "3′ processing and unsealed 5′-tail annealing · {} fragments · {:?} · homology {:?} bp",
            fragments.len(),
            policy.topology,
            policy.overlaps
        ),
        Operation::ProfiledOverlap {
            fragments,
            circular,
            options,
            ..
        } => format!(
            "Ordered {:?} overlap · {} fragments · {}",
            options.method,
            fragments.len(),
            if *circular { "circular" } else { "linear" }
        ),
        Operation::Overlap {
            fragments,
            circular,
            ..
        } => format!(
            "Ordered exact overlap · {} fragments · {}",
            fragments.len(),
            if *circular { "circular" } else { "linear" }
        ),
        Operation::Ligate {
            fragments,
            circular,
            allow_nicked,
            ..
        } => format!(
            "Ligation with bounded alternatives · {} fragments · {} · retained nicks {}",
            fragments.len(),
            if *circular { "circular" } else { "linear" },
            if *allow_nicked {
                "permitted"
            } else {
                "not permitted"
            }
        ),
        _ => "Other typed preparation; inspect the complete saved plan and protocol.".into(),
    }
}
pub(crate) struct Snapshot {
    pub(crate) plan: DesignPlan,
    pub(crate) target: Molecule,
    pub(crate) report: VerificationReport,
    pub(crate) plan_revision: String,
    pub(crate) target_revision: String,
    pub(crate) product_revision: Option<String>,
    pub(crate) error: Option<String>,
    pub(crate) physical: Option<DuplexSummary>,
    /// Bound leading documents for readonly inspection, never final/material authority.
    pub(crate) prefix_available: usize,
    pub(crate) prefix_error: Option<String>,
}
impl Snapshot {
    pub(crate) fn new(plan: DesignPlan, target: Molecule, report: VerificationReport) -> Self {
        let validation = (|| -> Result<(String, String), String> {
            collagene_core::assembly_junctions::source_preflight(&plan, &target)
                .map_err(|e| e.to_string())?;
            let plan_revision = verification::plan_revision(&plan);
            let target_revision = target.revision();
            if plan_revision.is_empty()
                || report.plan_revision != plan_revision
                || report.target_revision != target_revision
                || plan.target_revision != target_revision
            {
                return Err("Plan or independent-target revision does not match this reconstruction receipt.".into());
            }
            if report.engine_version != verification::ENGINE_VERSION {
                return Err("This reconstruction uses another scientific engine; reconstruct the current plan before inspecting its products.".into());
            }
            if report.products.len() > plan.steps.len() {
                return Err(
                    "The receipt contains more step products than the current plan.".into(),
                );
            }
            if report.status == Status::ConsistentUnderAssumptions
                && report.products.len() != plan.steps.len()
            {
                return Err(
                    "A consistent reconstruction receipt cannot omit step products.".into(),
                );
            }
            let products = report.products.iter().collect::<Vec<_>>();
            collagene_model::bounds::model_preflight(&products).map_err(|e| e.to_string())?;
            for molecule in &report.products {
                molecule.validate().map_err(|e| e.to_string())?;
                if molecule.metadata.get("plan_revision") != Some(&plan_revision)
                    || molecule.metadata.get("engine_version").map(String::as_str)
                        != Some(verification::ENGINE_VERSION)
                {
                    return Err(
                        "A step document belongs to a different plan or engine receipt.".into(),
                    );
                }
            }
            Ok((plan_revision, target_revision))
        })();
        let (plan_revision, target_revision, prefix_available, prefix_error) = match validation {
            Ok((p, t)) => (p, t, report.products.len(), None),
            Err(error) => (
                report.plan_revision.clone(),
                report.target_revision.clone(),
                0,
                Some(error),
            ),
        };
        let product_revision = if prefix_error.is_none() && prefix_available == plan.steps.len() {
            report.products.last().map(Molecule::revision)
        } else {
            None
        };
        let error=prefix_error.clone().or_else(||product_revision.is_none().then(||format!("Final product unavailable: {} of {} step products were reconstructed. Leading preparations are readonly diagnostic results, not the final assembly.",prefix_available,plan.steps.len())));
        let mut result = Self {
            plan,
            target,
            report,
            plan_revision,
            target_revision,
            product_revision,
            error,
            physical: None,
            prefix_available,
            prefix_error,
        };
        if result.product_revision.is_some() {
            result.physical = result.descriptive_physical_at(result.prefix_available - 1);
        }
        result
    }
    pub(crate) fn product_at(&self, index: usize) -> Result<&Molecule, String> {
        if let Some(error) = &self.prefix_error {
            return Err(error.clone());
        }
        if index >= self.prefix_available {
            return Err(
                "This step has no reconstructed product in the validated leading prefix.".into(),
            );
        }
        self.report
            .products
            .get(index)
            .ok_or_else(|| "The selected step product is unavailable.".into())
    }
    /// Descriptive metadata only; never restores tracked chemistry or authority.
    pub(crate) fn descriptive_physical_at(&self, index: usize) -> Option<DuplexSummary> {
        let product = self.product_at(index).ok()?;
        product
            .metadata
            .get("duplex_model")
            .filter(|text| text.len() <= 1024 * 1024)
            .and_then(|text| serde_json::from_str::<DuplexSummary>(text).ok())
            .filter(|summary| {
                summary.top_length == product.sequence.len()
                    && summary.topology == product.topology
                    && summary.nicks.len() <= 10000
                    && summary
                        .ends
                        .as_ref()
                        .is_none_or(|(a, b)| a.overhang.len() <= 10000 && b.overhang.len() <= 10000)
            })
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::{
        analysis::Primer, catalog_digest, tailing::TailingBase, verification::PlanInput,
    };
    use collagene_model::Topology;
    fn fixture(blocked: bool) -> (DesignPlan, Molecule, VerificationReport) {
        let source =
            Molecule::dna("Original parent", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
        let target = Molecule::dna(
            "Independent digest intent",
            "AATTCACTACCCCGGGGTAACG",
            Topology::Linear,
        )
        .unwrap();
        let catalog = catalog_digest::builtin_catalog().unwrap();
        let id = catalog
            .edition()
            .enzymes
            .iter()
            .find(|enzyme| enzyme.name == "EcoRI")
            .unwrap()
            .id
            .clone();
        let mut steps = vec![
            Operation::PcrWithMismatches {
                template: Source::Input("original".into()),
                forward: Primer::new("Mutation F", "ACTACCCC", "CCGAATTC").unwrap(),
                reverse: Primer::new("R", "GTTACCCC", "GGGAATTC").unwrap(),
                options: Default::default(),
            },
            Operation::CatalogDigest {
                template: Source::Step(0),
                selection: catalog_digest::select(&catalog, &[id]).unwrap(),
                fragment_index: 1,
            },
        ];
        if blocked {
            steps.push(Operation::Tail3Prime {
                template: Source::Step(1),
                base: TailingBase::A,
            });
        }
        let plan = DesignPlan {
            name: "Mismatch PCR and actual digest".into(),
            target_revision: target.revision(),
            inputs: vec![PlanInput {
                id: "original".into(),
                revision: source.revision(),
                molecule: source,
            }],
            steps,
        };
        let report = verification::verify(&plan, &target);
        (plan, target, report)
    }
    #[test]
    fn successfully_reconstructed_prefix_survives_unsupported_next_step_without_becoming_final() {
        let (plan, target, report) = fixture(true);
        assert_eq!(report.products.len(), 2, "{:?}", report.findings);
        assert_eq!(report.status, Status::Incomplete);
        let snapshot = Snapshot::new(plan, target, report);
        assert_eq!(snapshot.prefix_available, 2);
        assert!(snapshot.prefix_error.is_none());
        assert!(snapshot.product_revision.is_none() && snapshot.error.is_some());
        assert_eq!(
            snapshot.product_at(0).unwrap().sequence.as_str(),
            "CCGAATTCACTACCCCGGGGTAACGAATTCCC"
        );
        assert_eq!(
            snapshot.product_at(1).unwrap().sequence.as_str(),
            "AATTCACTACCCCGGGGTAACG"
        );
        assert!(snapshot.product_at(2).is_err());
        assert!(snapshot.descriptive_physical_at(1).is_some());
        assert!(snapshot.physical.is_none());
    }
    #[test]
    fn complete_native_result_retains_existing_final_availability() {
        let (plan, target, report) = fixture(false);
        assert_eq!(
            report.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            report.findings
        );
        let snapshot = Snapshot::new(plan, target, report);
        assert_eq!(snapshot.prefix_available, 2);
        assert_eq!(
            snapshot.product_revision.as_deref(),
            Some(snapshot.product_at(1).unwrap().revision().as_str())
        );
        assert!(snapshot.error.is_none() && snapshot.physical.is_some());
    }
    #[test]
    fn stale_shared_context_or_malformed_counts_withhold_every_prefix() {
        let (plan, target, report) = fixture(true);
        let mut changed = plan.clone();
        changed.inputs[0].id = "different-role".into();
        let s = Snapshot::new(changed, target.clone(), report.clone());
        assert_eq!(s.prefix_available, 0);
        assert!(s.product_at(0).is_err());
        let mut changed = target.clone();
        changed
            .metadata
            .insert("changed".into(), "target identity".into());
        assert_eq!(
            Snapshot::new(plan.clone(), changed, report.clone()).prefix_available,
            0
        );
        let mut changed = report.clone();
        changed.plan_revision = "0".repeat(64);
        assert_eq!(
            Snapshot::new(plan.clone(), target.clone(), changed).prefix_available,
            0
        );
        let mut changed = report.clone();
        changed.engine_version = "historical".into();
        assert_eq!(
            Snapshot::new(plan.clone(), target.clone(), changed).prefix_available,
            0
        );
        let mut changed = report.clone();
        changed.products[0]
            .metadata
            .insert("plan_revision".into(), "0".repeat(64));
        assert_eq!(
            Snapshot::new(plan.clone(), target.clone(), changed).prefix_available,
            0
        );
        let mut changed = report.clone();
        changed.products[1]
            .metadata
            .insert("engine_version".into(), "historical".into());
        assert_eq!(
            Snapshot::new(plan.clone(), target.clone(), changed).prefix_available,
            0
        );
        let mut changed = report.clone();
        changed.status = Status::ConsistentUnderAssumptions;
        assert_eq!(
            Snapshot::new(plan.clone(), target.clone(), changed).prefix_available,
            0
        );
        let mut changed = report;
        changed.products.resize(4, changed.products[0].clone());
        assert_eq!(Snapshot::new(plan, target, changed).prefix_available, 0);
    }
    #[test]
    fn descriptive_metadata_never_changes_context_or_final_availability() {
        let (plan, target, mut report) = fixture(true);
        report.products[0]
            .metadata
            .insert("duplex_model".into(), "{}".into());
        let snapshot = Snapshot::new(plan, target, report);
        assert_eq!(snapshot.prefix_available, 2);
        assert!(snapshot.descriptive_physical_at(0).is_none());
        assert!(snapshot.product_revision.is_none());
        assert!(operation_label(&snapshot.plan.steps[0]).contains("Intentional-mismatch PCR"));
        assert!(operation_label(&snapshot.plan.steps[2]).contains("3′ tails"));
    }
}
