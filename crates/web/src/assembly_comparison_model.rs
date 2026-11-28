//! Diagnostic comparison of the independently intended reference and actual final.
//! No history, document, target replacement or physical material is constructed.
use crate::assembly_product_model::Snapshot;
use collagene_core::{
    pairwise_review::{self, Options, Report},
    verification,
};
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Capture {
    pub(crate) context_revision: String,
    pub(crate) plan_revision: String,
    pub(crate) target_revision: String,
    pub(crate) product_revision: String,
    pub(crate) pairwise_revision: String,
    pub(crate) options: Options,
}
pub(crate) fn capture(snapshot: &Snapshot, options: Options) -> Result<Capture, String> {
    if snapshot.prefix_error.is_some()
        || snapshot.product_revision.is_none()
        || snapshot.report.products.len() != snapshot.plan.steps.len()
        || snapshot.plan.steps.is_empty()
    {
        return Err("A bound actual final product is required; a leading preparation prefix cannot be compared as final.".into());
    }
    collagene_core::assembly_junctions::source_preflight(&snapshot.plan, &snapshot.target)
        .map_err(|e| e.to_string())?;
    let plan_revision = verification::plan_revision(&snapshot.plan);
    let target_revision = snapshot.target.revision();
    if plan_revision != snapshot.plan_revision
        || plan_revision != snapshot.report.plan_revision
        || target_revision != snapshot.target_revision
        || target_revision != snapshot.plan.target_revision
        || target_revision != snapshot.report.target_revision
        || snapshot.report.engine_version != verification::ENGINE_VERSION
    {
        return Err(
            "The reconstruction plan or independent target changed; obtain a new snapshot.".into(),
        );
    }
    let product = snapshot.product_at(snapshot.plan.steps.len() - 1)?;
    let product_revision = product.revision();
    if snapshot.product_revision.as_ref() != Some(&product_revision)
        || product.metadata.get("plan_revision") != Some(&plan_revision)
        || product.metadata.get("engine_version").map(String::as_str)
            != Some(verification::ENGINE_VERSION)
    {
        return Err("The actual final-product receipt changed; reconstruct it again.".into());
    }
    let pairwise_revision = pairwise_review::revision(&snapshot.target, product, &options)
        .map_err(|e| e.to_string())?;
    let bytes = serde_json::to_vec(&(
        "assembly-final-comparison-1",
        &plan_revision,
        &target_revision,
        &product_revision,
        &pairwise_revision,
        options,
    ))
    .map_err(|e| e.to_string())?;
    let context_revision = blake3::hash(&bytes).to_hex().to_string();
    Ok(Capture {
        context_revision,
        plan_revision,
        target_revision,
        product_revision,
        pairwise_revision,
        options,
    })
}
pub(crate) fn validate(
    snapshot: &Snapshot,
    bound: &Capture,
    report: &Report,
) -> Result<(), String> {
    if &capture(snapshot, bound.options)? != bound {
        return Err("The plan, target, actual final or comparison policy changed.".into());
    }
    if report.options != bound.options
        || report.reference_revision != bound.target_revision
        || report.query_revision != bound.product_revision
    {
        return Err(
            "Comparison sources or options differ from the captured final reconstruction.".into(),
        );
    }
    report
        .validate_binding(&bound.pairwise_revision)
        .map_err(|e| e.to_string())?;
    let actual = snapshot.product_at(snapshot.plan.steps.len() - 1)?;
    report
        .validate_sources(&snapshot.target, actual)
        .map_err(|e| e.to_string())
}
/// Native alignment coordinates map back to original reference/query coordinates,
/// preserving reverse orientation, circular origin crossings and insertion carets.
pub(crate) fn selected_difference(
    snapshot: &Snapshot,
    bound: &Capture,
    report: &Report,
    index: usize,
) -> Result<(Vec<(usize, usize)>, (usize, usize)), String> {
    validate(snapshot, bound, report)?;
    let difference = report
        .comparison
        .alignment
        .differences
        .get(index)
        .ok_or("The selected difference is unavailable.")?;
    let reference = report
        .map_reference(difference.reference_start, difference.reference_end)
        .map_err(|e| e.to_string())?;
    let query = report
        .map_query(difference.query_start, difference.query_end)
        .map_err(|e| e.to_string())?;
    Ok((reference, query))
}
#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::{
        alignment::{ComparisonOptions, DifferenceKind},
        analysis::Primer,
        verification::{DesignPlan, Operation, PlanInput, Source, Status},
    };
    use collagene_model::{Molecule, Topology};
    fn snapshot(target: &str, circular: bool) -> Snapshot {
        let source =
            Molecule::dna("Original source", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
        let target = Molecule::dna(
            "Independent intended product",
            target,
            if circular {
                Topology::Circular
            } else {
                Topology::Linear
            },
        )
        .unwrap();
        let plan = DesignPlan {
            name: "Exact PCR diagnostic".into(),
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
        let report = verification::verify(&plan, &target);
        Snapshot::new(plan, target, report)
    }
    fn compare(snapshot: &Snapshot, options: Options) -> (Capture, Report) {
        let bound = capture(snapshot, options).unwrap();
        let report = pairwise_review::compare(
            &snapshot.target,
            snapshot.report.products.last().unwrap(),
            &options,
        )
        .unwrap();
        validate(snapshot, &bound, &report).unwrap();
        (bound, report)
    }
    #[test]
    fn invalid_actual_final_can_be_diagnosed_against_unchanged_intended_snp_and_indel() {
        let options = Options::Stored(Default::default());
        for (target, kind) in [
            ("ACTACCCCGGGGTAAC", DifferenceKind::Substitution),
            ("ACGATCCCCGGGGTAAC", DifferenceKind::Deletion),
        ] {
            let snapshot = snapshot(target, false);
            assert_eq!(snapshot.report.status, Status::Invalid);
            let original = snapshot.target.revision();
            let (bound, report) = compare(&snapshot, options);
            assert_eq!(report.reference_revision, snapshot.target_revision);
            assert_eq!(
                report.query_revision,
                snapshot.product_revision.clone().unwrap()
            );
            let index = report
                .comparison
                .alignment
                .differences
                .iter()
                .position(|difference| difference.kind == kind)
                .unwrap();
            let (reference, query) =
                selected_difference(&snapshot, &bound, &report, index).unwrap();
            let difference = &report.comparison.alignment.differences[index];
            assert_eq!(
                reference,
                vec![(difference.reference_start, difference.reference_end)]
            );
            assert_eq!(query, (difference.query_start, difference.query_end));
            assert_eq!(snapshot.target.revision(), original);
        }
    }
    #[test]
    fn whole_comparison_rotation_and_reverse_are_diagnostic_not_verifier_reinterpretation() {
        // The actual PCR is linear. Circular target topology and reversed target DNA
        // remain verifier mismatches even when whole-reference alignment finds identity.
        for target in ["CCCCGGGGTAACACGA", "GTTACCCCGGGGTCGT"] {
            let snapshot = snapshot(target, true);
            assert_eq!(snapshot.report.status, Status::Invalid);
            let (_, report) = compare(&snapshot, Options::Whole(ComparisonOptions::default()));
            assert!(report.comparison.alignment.differences.is_empty());
            assert!(report.comparison.reference_origin != 0 || report.comparison.query_reversed);
            assert_eq!(snapshot.report.status, Status::Invalid);
        }
    }
    #[test]
    fn selected_differences_map_rotated_and_reverse_searches_back_to_original_bases() {
        for (target, reference_position) in [("CCCCGGGGTAACACTA", 14), ("GTTACCCCGGGGTAGT", 13)] {
            let snapshot = snapshot(target, true);
            let (bound, report) = compare(&snapshot, Options::Whole(Default::default()));
            let index = report
                .comparison
                .alignment
                .differences
                .iter()
                .position(|difference| difference.kind == DifferenceKind::Substitution)
                .unwrap();
            let (reference, query) =
                selected_difference(&snapshot, &bound, &report, index).unwrap();
            assert_eq!(
                reference,
                vec![(reference_position, reference_position + 1)]
            );
            assert_eq!(query, (2, 3));
        }
    }
    #[test]
    fn stale_snapshot_options_or_reply_binding_cannot_retarget_comparison() {
        let snapshot = snapshot("ACTACCCCGGGGTAAC", false);
        let (bound, report) = compare(&snapshot, Options::Stored(Default::default()));
        let mut changed = bound.clone();
        changed.plan_revision = "0".repeat(64);
        assert!(validate(&snapshot, &changed, &report).is_err());
        let mut changed = bound.clone();
        changed.options = Options::Whole(Default::default());
        assert!(validate(&snapshot, &changed, &report).is_err());
        let mut changed = report.clone();
        changed.query_revision = changed.reference_revision.clone();
        assert!(validate(&snapshot, &bound, &changed).is_err());
        let mut changed = report;
        changed.options = Options::Whole(Default::default());
        assert!(validate(&snapshot, &bound, &changed).is_err());
        let mut stale = snapshot;
        stale.plan.name.push_str(" edited");
        assert!(capture(&stale, bound.options).is_err());
    }
    #[test]
    fn partial_prefix_cannot_be_compared_as_a_final_and_postcapture_product_changes_fail() {
        let mut snapshot = snapshot("ACTACCCCGGGGTAAC", false);
        let (bound, report) = compare(&snapshot, Options::Stored(Default::default()));
        snapshot.report.products[0].name.push_str(" changed");
        assert!(validate(&snapshot, &bound, &report).is_err());
        let initial = self::snapshot("ACTACCCCGGGGTAAC", false);
        let mut plan = initial.plan;
        plan.steps.push(Operation::Digest {
            template: Source::Step(0),
            enzymes: vec!["EcoRI".into()],
            fragment_index: 999,
        });
        let report = verification::verify(&plan, &initial.target);
        let partial = Snapshot::new(plan, initial.target, report);
        assert_eq!(partial.prefix_available, 1);
        assert!(partial.product_revision.is_none());
        assert!(capture(&partial, Options::Stored(Default::default())).is_err());
        let mut snapshot = self::snapshot("ACTACCCCGGGGTAAC", false);
        snapshot.target.name.push_str(" changed");
        assert!(capture(&snapshot, Options::Stored(Default::default())).is_err());
    }
}
