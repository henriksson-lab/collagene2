//! Actual scientific work owned by one plan or one causal-QC replay job.
//! This state is borrowed in process, never deserialized from an assessment.
use crate::{chemistry::tracked_duplex::CleavageBudget, sequence::analysis::ExactSearchBudget};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum ReplayScope {
    /// Existing verification behavior: nested retained producers keep their own
    /// original native scientific scopes.
    PlanLocal,
    /// A causal-QC request shares actual exact-PCR and cleavage charges across every
    /// retained producer it replays, including nested saved material.
    QcJob,
}

pub(crate) struct ReplayBudget<'a> {
    pub(crate) pcr: ExactSearchBudget,
    pub(crate) cleavage: &'a mut CleavageBudget,
    scope: ReplayScope,
}
impl<'a> ReplayBudget<'a> {
    pub fn plan_local(cleavage: &'a mut CleavageBudget) -> Self {
        Self {
            pcr: ExactSearchBudget::default(),
            cleavage,
            scope: ReplayScope::PlanLocal,
        }
    }
    pub fn qc_job(cleavage: &'a mut CleavageBudget) -> Self {
        Self {
            pcr: ExactSearchBudget::default(),
            cleavage,
            scope: ReplayScope::QcJob,
        }
    }
    pub fn is_qc_job(&self) -> bool {
        self.scope == ReplayScope::QcJob
    }
    pub fn exhausted(&self) -> bool {
        self.pcr.remaining() == 0 || self.cleavage.remaining() == 0
    }
    pub fn pcr_used(&self) -> usize {
        self.pcr.used()
    }
    pub fn cleavage_used(&self) -> usize {
        self.cleavage.used()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        chemistry::Strand,
        planning::plan::{DesignPlan, Operation, PlanInput, Source},
        sequence::analysis::Primer,
        verification::engine,
    };
    use collagene_model::{Molecule, Topology};

    fn fixture() -> (DesignPlan, Molecule) {
        let source = Molecule::dna(
            "literal template",
            "AGTCCGATGCTACGTAGGACTTCAGACCTGTA",
            Topology::Linear,
        )
        .unwrap();
        let target = Molecule::dna(
            "literal expected PCR",
            "GGAGTCCGATGCTACGTAGGACTTCAGACCTGTACC",
            Topology::Linear,
        )
        .unwrap();
        let plan = DesignPlan {
            name: "exact PCR budget fixture".into(),
            target_revision: target.revision(),
            inputs: vec![PlanInput {
                id: "template".into(),
                revision: source.revision(),
                molecule: source,
            }],
            steps: vec![Operation::Pcr {
                template: Source::Input("template".into()),
                forward: Primer::new("forward", "AGTCCGAT", "GG").unwrap(),
                reverse: Primer::new("reverse", "TACAGGTC", "GG").unwrap(),
                template_strand: Strand::Top,
            }],
        };
        (plan, target)
    }
    #[test]
    fn separate_native_replays_share_actual_qc_pcr_work_even_without_trace_capture() {
        let (plan, target) = fixture();
        let mut cleavage = CleavageBudget::new(100_000_000).unwrap();
        let mut baseline = ReplayBudget::plan_local(&mut cleavage);
        let (report, _) = engine::verify_with_tracked_scoped(&plan, &target, &mut baseline);
        assert_eq!(report.product(0).unwrap().sequence, target.sequence);
        let first_work = baseline.pcr_used();
        assert!(first_work > 0);

        let mut cleavage = CleavageBudget::new(100_000_000).unwrap();
        let mut qc = ReplayBudget::qc_job(&mut cleavage);
        // Lower the private test ledger instead of constructing an enormous
        // sequence. Each actual kernel charge still follows the production path.
        let limit = first_work + first_work / 2;
        qc.pcr = ExactSearchBudget::new(limit).unwrap();
        let (first, _) = engine::verify_with_tracked_scoped(&plan, &target, &mut qc);
        assert!(first.product(0).is_some());
        assert_eq!(qc.pcr_used(), first_work);
        let (second, _) = engine::verify_with_tracked_scoped(&plan, &target, &mut qc);
        assert!(second.product(0).is_none());
        assert_eq!(qc.pcr_used(), limit);
        assert!(qc.exhausted());
        let (third, _) = engine::verify_with_tracked_scoped(&plan, &target, &mut qc);
        assert!(third.product(0).is_none());
        assert_eq!(qc.pcr_used(), limit);

        // A separate ordinary verification retains its own existing plan scope.
        assert!(engine::verify(&plan, &target).product(0).is_some());
    }
}
