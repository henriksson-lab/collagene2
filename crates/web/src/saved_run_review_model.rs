//! Exact saved-run review scope; current biological heads are never source selectors.
use collagene_core::{
    run_review::PreparedReview,
    verification::{self, VerificationReport},
};
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Scope {
    pub epoch: u64,
    pub workspace_fingerprint: String,
    pub run_id: String,
    pub run_attachment_hash: String,
}
impl Scope {
    pub fn permit(&self, current: &Self) -> Result<(), String> {
        if self != current {
            return Err(
                "The workspace or exact saved-run attachment changed. Prepare this review again."
                    .into(),
            );
        }
        Ok(())
    }
    pub fn validate_prepared(&self, p: &PreparedReview) -> Result<(), String> {
        if p.source_workspace_fingerprint != self.workspace_fingerprint
            || p.run_id != self.run_id
            || p.run_attachment_hash != self.run_attachment_hash
        {
            return Err("The prepared review belongs to another workspace or saved run.".into());
        }
        Ok(())
    }
}
/// Bind an editing handoff to the unchanged exact native review, not current heads.
pub fn editing_candidate(
    expected: &Scope,
    current: &Scope,
    prepared: &PreparedReview,
    workspace: &collagene_core::workspace::Workspace,
) -> Result<collagene_core::transfer::Candidate, String> {
    expected.permit(current)?;
    expected.validate_prepared(prepared)?;
    prepared
        .validate_for_workspace(workspace)
        .map_err(|e| e.to_string())?;
    Ok(prepared.candidate.clone())
}

pub fn validate_report(p: &PreparedReview, r: &VerificationReport) -> Result<(), String> {
    if r.engine_version != verification::ENGINE_VERSION
        || r.plan_revision != p.candidate.plan_revision
        || r.target_revision != p.target_event.revision
        || r.target_revision != p.target.revision()
        || r.products.len() > p.candidate.plan.steps.len()
    {
        return Err(
            "Fresh verification does not match the exact saved plan and retained intended target."
                .into(),
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::{
        Molecule, Topology,
        analysis::Primer,
        run_lineage::{self, EventBinding, InputEventBinding},
        run_review,
        verification::{DesignPlan, Operation, PlanInput, Source, Status},
        workspace::{RevisionEvent, Workspace},
    };
    fn fixture() -> (Workspace, PreparedReview, VerificationReport) {
        let source = Molecule::dna("Stock", "AGTCCGATGCTACGTA", Topology::Linear).unwrap();
        let target =
            Molecule::dna("Independent intent", "AGTCCGATGCTACGTA", Topology::Linear).unwrap();
        let w =
            Workspace::from_documents("Review", vec![source.clone(), target.clone()], 1).unwrap();
        let pin = |i: usize| {
            let d = &w.documents[i];
            EventBinding {
                document_id: d.id.clone(),
                event_id: d.head_event.clone().unwrap(),
                revision: d.head.clone(),
            }
        };
        let plan = DesignPlan {
            name: "PCR".into(),
            target_revision: target.revision(),
            inputs: vec![PlanInput {
                id: "stock".into(),
                revision: source.revision(),
                molecule: source,
            }],
            steps: vec![Operation::Pcr {
                template_strand: collagene_core::chemistry::Strand::Top,
                template: Source::Input("stock".into()),
                forward: Primer::new("F", "AGTCC", "").unwrap(),
                reverse: Primer::new("R", "TACGT", "").unwrap(),
            }],
        };
        let saved = run_lineage::record_run(
            &w,
            &plan,
            pin(1),
            vec![InputEventBinding {
                input_id: "stock".into(),
                event: pin(0),
            }],
            false,
        )
        .unwrap();
        let p = run_review::prepare(&saved.workspace, &saved.receipt.run_id).unwrap();
        let report = verification::verify(&p.candidate.plan, &p.target);
        (saved.workspace, p, report)
    }
    fn bound(w: &Workspace, p: &PreparedReview, epoch: u64) -> Scope {
        Scope {
            epoch,
            workspace_fingerprint: run_lineage::workspace_fingerprint(w).unwrap(),
            run_id: p.run_id.clone(),
            run_attachment_hash: p.run_attachment_hash.clone(),
        }
    }
    #[test]
    fn fresh_report_must_bind_both_saved_plan_and_independent_target() {
        let (_, p, r) = fixture();
        validate_report(&p, &r).unwrap();
        let other =
            Molecule::dna("Another independent target", "TTTTCCCC", Topology::Linear).unwrap();
        let wrong = verification::verify(&p.candidate.plan, &other);
        assert!(validate_report(&p, &wrong).is_err());
        let mut changed = p.candidate.plan.clone();
        changed.name = "Edited declaration".into();
        assert!(validate_report(&p, &verification::verify(&changed, &p.target)).is_err());
    }
    #[test]
    fn equal_dna_new_event_invalidates_the_export_scope() {
        let (mut w, p, _) = fixture();
        let original = bound(&w, &p, 7);
        let d = &mut w.documents[0];
        let e = RevisionEvent::new(
            d.id.clone(),
            d.head.clone(),
            vec![d.head_event.clone().unwrap()],
            "Reopened equal DNA under a new event",
        )
        .unwrap();
        d.head_event = Some(e.id.clone());
        w.history.push(e);
        w.validate().unwrap();
        assert!(original.permit(&bound(&w, &p, 7)).is_err());
        let mut changed = original.clone();
        changed.epoch += 1;
        assert!(original.permit(&changed).is_err());
    }
    #[test]
    fn diagnostic_status_is_current_evidence_not_historical_acceptance() {
        let (_, mut p, mut r) = fixture();
        p.recorded_assessment.status = Status::ConsistentUnderAssumptions;
        r.status = Status::Incomplete;
        validate_report(&p, &r).unwrap();
        r.engine_version = "historical-engine/science-1".into();
        assert!(validate_report(&p, &r).is_err());
        let mut scope = Scope {
            epoch: 1,
            workspace_fingerprint: p.source_workspace_fingerprint.clone(),
            run_id: p.run_id.clone(),
            run_attachment_hash: p.run_attachment_hash.clone(),
        };
        scope.validate_prepared(&p).unwrap();
        scope.run_attachment_hash = "00".repeat(32);
        assert!(scope.validate_prepared(&p).is_err());
    }
    #[test]
    fn editing_handoff_keeps_retained_events_and_rejects_changed_scope() {
        let (workspace, prepared, _) = fixture();
        let scope = bound(&workspace, &prepared, 12);
        let candidate = editing_candidate(&scope, &scope, &prepared, &workspace).unwrap();
        assert_eq!(candidate.plan_revision, prepared.candidate.plan_revision);
        assert_eq!(
            candidate.plan.target_revision,
            prepared.target_event.revision
        );
        assert_eq!(
            run_lineage::workspace_fingerprint(&workspace).unwrap(),
            scope.workspace_fingerprint
        );
        let mut later = scope.clone();
        later.epoch += 1;
        assert!(editing_candidate(&scope, &later, &prepared, &workspace).is_err());
        let mut forged = prepared.clone();
        forged.inputs[0].event.event_id = "00".repeat(32);
        assert!(editing_candidate(&scope, &scope, &forged, &workspace).is_err());
    }
}
