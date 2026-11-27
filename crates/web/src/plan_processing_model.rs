//! Raw physical-processing declarations shared by the candidate editor and
//! native tests. Parsing never grants physical substrate or terminal authority.
use collagene_core::{
    three_prime_resection::Cuts,
    verification::{DesignPlan, Operation, Source},
};

#[derive(Clone, Default, PartialEq, Eq)]
pub(crate) struct ResectionDraft {
    pub left: String,
    pub right: String,
}
impl ResectionDraft {
    pub fn capture(cuts: &Cuts) -> Self {
        Self {
            left: cuts.left.to_string(),
            right: cuts.right.to_string(),
        }
    }
    pub fn cuts(&self) -> Result<Cuts, String> {
        fn count(raw: &str) -> Result<usize, String> {
            let text = raw.trim();
            if raw.len() > 20 || text.is_empty() || !text.bytes().all(|b| b.is_ascii_digit()) {
                return Err("Enter an explicit nonnegative whole number for each 3′ cut, including zero for an end left unchanged.".into());
            }
            text.parse::<usize>()
                .ok()
                .filter(|n| *n <= collagene_core::three_prime_resection::MAX_SOURCE_BASES / 2)
                .ok_or_else(|| "A 3′ cut exceeds the native limit of 1,000,000 bases.".into())
        }
        Ok(Cuts {
            left: count(&self.left)?,
            right: count(&self.right)?,
        })
    }
    pub fn apply(&self, operation: &Operation) -> Result<Operation, String> {
        let Operation::ResectThreePrime { source, .. } = operation else {
            return Err(
                "3′ cut settings belong only to an explicit 3′ source-resection operation.".into(),
            );
        };
        Ok(Operation::ResectThreePrime {
            source: source.clone(),
            cuts: self.cuts()?,
        })
    }
}

/// Static documentary ancestry only. This does not grant material authority to
/// unflagged steps: fresh native reconstruction still decides reusability.
pub(crate) fn reference_only_steps(plan: &DesignPlan) -> Vec<bool> {
    let mut flags = Vec::with_capacity(plan.steps.len());
    for operation in &plan.steps {
        let direct = matches!(
            operation,
            Operation::Topo { .. }
                | Operation::CheckedTopo { .. }
                | Operation::RecombineReference { .. }
                | Operation::MultiSiteRecombineReference { .. }
                | Operation::Overlap { .. }
                | Operation::ProfiledOverlap { .. }
                | Operation::HomologyAnnealing { .. }
        );
        let inherited = collagene_core::verification::operation_sources(operation)
            .into_iter()
            .any(|source| matches!(source,Source::Step(i) if flags.get(*i)==Some(&true)));
        flags.push(direct || inherited);
    }
    flags
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::{
        Molecule, Topology,
        analysis::Primer,
        chemistry::ThreePrime,
        verification::{self, PhysicalProcessingEvidence, PlanInput, Status},
    };
    const STOCK: &str = "AGTCCGATGCTACGTAGGACTTCA";
    const TRIMMED: &str = "AGTCCGATGCTACGTAGGA";
    fn fixture() -> (DesignPlan, Molecule) {
        let stock = Molecule::dna("Original stock", STOCK, Topology::Linear).unwrap();
        let target = Molecule::dna(
            "Independent shortened top reference",
            TRIMMED,
            Topology::Linear,
        )
        .unwrap();
        (
            DesignPlan {
                name: "PCR then explicit3′ cuts".into(),
                target_revision: target.revision(),
                inputs: vec![PlanInput {
                    id: "stock".into(),
                    revision: stock.revision(),
                    molecule: stock,
                }],
                steps: vec![Operation::Pcr {
                    template_strand: collagene_core::chemistry::Strand::Top,
                    template: Source::Input("stock".into()),
                    forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
                    reverse: Primer::new("R", "TGAAGTCC", "").unwrap(),
                }],
            },
            target,
        )
    }
    #[test]
    fn explicit_cut_edit_preserves_dependency_and_replays_real_pcr_strands() {
        let (mut plan, target) = fixture();
        let raw = ResectionDraft {
            left: "3".into(),
            right: "5".into(),
        };
        let original = Operation::ResectThreePrime {
            source: Source::Step(0),
            cuts: Cuts { left: 0, right: 0 },
        };
        let edited = raw.apply(&original).unwrap();
        assert!(
            matches!(&edited,Operation::ResectThreePrime{source:Source::Step(0),cuts} if cuts.left==3&&cuts.right==5)
        );
        plan.steps.push(edited);
        let report = verification::verify(&plan, &target);
        assert_eq!(
            report.status,
            Status::ConsistentUnderAssumptions,
            "{:?}",
            report.findings
        );
        assert_eq!(report.products.last().unwrap().sequence.as_str(), TRIMMED);
        let PhysicalProcessingEvidence::ThreePrimeResection(geometry) =
            &report.physical_processing_reports.last().unwrap().evidence
        else {
            panic!("Actual native resection evidence");
        };
        assert_eq!(geometry.product.top.as_str(), TRIMMED);
        assert_eq!(geometry.product.bottom.as_str(), "TGAAGTCCTACGTAGCATCGG");
        assert_eq!(geometry.product.bottom_start, 3);
        assert_eq!(geometry.product.left.three_prime, ThreePrime::Hydroxyl);
        assert_eq!(geometry.product.right.three_prime, ThreePrime::Hydroxyl);
        assert!(geometry.product.nicks.is_empty());
        // The same letters without their preceding PCR do not establish3′OH.
        plan.steps = vec![
            raw.apply(&Operation::ResectThreePrime {
                source: Source::Input("stock".into()),
                cuts: Cuts { left: 0, right: 0 },
            })
            .unwrap(),
        ];
        let unknown = verification::verify(&plan, &target);
        assert_eq!(unknown.status, Status::Incomplete);
        assert!(unknown.physical_processing_reports.is_empty());
        assert!(unknown.products.is_empty());
    }
    #[test]
    fn raw_cut_fields_preserve_zero_and_refuse_unfinished_or_unbounded_values() {
        let mut raw = ResectionDraft::default();
        assert!(raw.cuts().is_err());
        for bad in [
            "",
            "partial",
            "-1",
            "+1",
            "1.5",
            "1000001",
            "999999999999999999999",
        ] {
            raw.left = bad.into();
            raw.right = "0".into();
            assert!(raw.cuts().is_err());
        }
        raw = ResectionDraft::capture(&Cuts {
            left: 0,
            right: 1_000_000,
        });
        assert_eq!(
            raw.cuts().unwrap(),
            Cuts {
                left: 0,
                right: 1_000_000
            }
        );
        assert!(
            raw.apply(&Operation::Pcr {
                template_strand: collagene_core::chemistry::Strand::Top,
                template: Source::Input("stock".into()),
                forward: Primer::new("F", "ACGT", "").unwrap(),
                reverse: Primer::new("R", "TGCA", "").unwrap()
            })
            .is_err()
        );
    }
    #[test]
    fn annealed_reference_ancestry_propagates_to_downstream_preparations_only() {
        let (mut plan, _) = fixture();
        plan.steps.push(Operation::HomologyAnnealing {
            fragments: vec![verification::OrientedSource {
                source: Source::Step(0),
                reverse: false,
            }],
            options: Default::default(),
            policy: collagene_core::homology_annealing::Options {
                planned_rescue: None,
                overlaps: vec![4],
                topology: Topology::Circular,
                cuts: vec![Cuts { left: 4, right: 4 }],
                min_paired_bases: 4,
                max_comparisons: 100,
                max_candidates: 10,
            },
        });
        plan.steps.push(Operation::ResectThreePrime {
            source: Source::Step(1),
            cuts: Cuts { left: 1, right: 1 },
        });
        plan.steps.push(Operation::ResectThreePrime {
            source: Source::Step(0),
            cuts: Cuts { left: 1, right: 1 },
        });
        assert_eq!(reference_only_steps(&plan), vec![false, true, true, false]);
    }
}
