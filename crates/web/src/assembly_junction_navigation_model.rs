//! Exact navigation in the actual reconstructed final product. No history,
//! reference target, physical state or verification result is created here.
use crate::assembly_product_model::Snapshot;
use collagene_core::{
    assembly_junctions::{self, Options, Report},
    verification,
};

const ENGINE_VERSION: &str = "assembly-junction-navigation-1";

#[derive(Clone, Debug)]
pub(crate) struct Capture {
    pub(crate) context_revision: String,
    pub(crate) plan_revision: String,
    pub(crate) target_revision: String,
    pub(crate) product_revision: String,
    pub(crate) inspection_revision: String,
    pub(crate) options: Options,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Kind {
    Window,
    Homology,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Selection {
    pub(crate) context_revision: String,
    pub(crate) product_revision: String,
    pub(crate) junction_index: usize,
    pub(crate) kind: Kind,
    /// Zero-based half-open ranges in the original ACTUAL final product, in
    /// native biological window order. A circular origin crossing stays split.
    pub(crate) parts: Vec<(usize, usize)>,
}
fn options_equal(a: &Options, b: &Options) -> bool {
    a.flank_bases == b.flank_bases
        && a.genetic_code == b.genetic_code
        && a.frames == b.frames
        && a.max_junctions == b.max_junctions
        && a.include_cds == b.include_cds
}

/// Capture an already reconstructed final product. A partial preparation prefix
/// is never reinterpreted as final, including when its letters happen to agree.
/// The host must additionally guard workspace epoch and current UI ownership.
pub(crate) fn capture(snapshot: &Snapshot, options: Options) -> Result<Capture, String> {
    if snapshot.prefix_error.is_some()
        || snapshot.product_revision.is_none()
        || snapshot.plan.steps.is_empty()
        || snapshot.report.products.len() != snapshot.plan.steps.len()
        || snapshot.prefix_available != snapshot.plan.steps.len()
    {
        return Err("Junction navigation requires the actual final product, not an intermediate preparation.".into());
    }
    let inspection_revision =
        assembly_junctions::inspection_revision(&snapshot.plan, &snapshot.target, &options)
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
            "The original plan or independent target changed; reconstruct the current plan.".into(),
        );
    }
    let product = snapshot.product_at(snapshot.plan.steps.len() - 1)?;
    // Recheck borrowed bounds before hashing the full current molecule. Snapshot
    // validation is a construction receipt, not permission to trust later edits.
    collagene_core::model::bounds::model_preflight(&[product]).map_err(|e| e.to_string())?;
    product.validate().map_err(|e| e.to_string())?;
    if product.sequence.alphabet() != collagene_core::Alphabet::Dna {
        return Err("Junction navigation requires an actual reconstructed DNA product.".into());
    }
    let product_revision = product.revision();
    if snapshot.product_revision.as_ref() != Some(&product_revision)
        || product.metadata.get("plan_revision") != Some(&plan_revision)
        || product.metadata.get("engine_version").map(String::as_str)
            != Some(verification::ENGINE_VERSION)
    {
        return Err(
            "The full actual-product receipt changed; inspect a fresh final reconstruction.".into(),
        );
    }
    let bytes = serde_json::to_vec(&(
        ENGINE_VERSION,
        &plan_revision,
        &target_revision,
        &product_revision,
        &inspection_revision,
    ))
    .map_err(|e| e.to_string())?;
    Ok(Capture {
        context_revision: blake3::hash(&bytes).to_hex().to_string(),
        plan_revision,
        target_revision,
        product_revision,
        inspection_revision,
        options,
    })
}

/// Readonly binding/shape checks of the worker receipt against the captured
/// snapshots. This does not rerun molecular reconstruction or upgrade its verdict.
pub(crate) fn validate(
    snapshot: &Snapshot,
    bound: &Capture,
    report: &Report,
) -> Result<(), String> {
    let current = capture(snapshot, bound.options.clone())?;
    if current.context_revision != bound.context_revision
        || current.plan_revision != bound.plan_revision
        || current.target_revision != bound.target_revision
        || current.product_revision != bound.product_revision
        || current.inspection_revision != bound.inspection_revision
        || !options_equal(&current.options, &bound.options)
    {
        return Err("The captured plan, target, product or junction settings changed.".into());
    }
    report
        .validate_binding(&bound.inspection_revision)
        .map_err(|e| e.to_string())?;
    if report.plan_revision != bound.plan_revision
        || report.target_revision != bound.target_revision
        || !options_equal(&report.options, &bound.options)
        || report.verification.engine_version != verification::ENGINE_VERSION
        || report.verification.products.len() != snapshot.plan.steps.len()
        || report.product_revision.as_ref() != Some(&bound.product_revision)
    {
        return Err(
            "Junction inspection does not belong to this actual final-product snapshot.".into(),
        );
    }
    let actual = snapshot.product_at(snapshot.plan.steps.len() - 1)?;
    let inspected = report
        .verification
        .products
        .last()
        .ok_or("Junction inspection has no actual final product.")?;
    if inspected != actual
        || inspected.revision() != bound.product_revision
        || report.product_length != Some(actual.sequence.len())
        || report.product_topology != Some(actual.topology)
    {
        return Err(
            "The complete inspected product differs from the current final product.".into(),
        );
    }
    Ok(())
}

/// Return only native product-coordinate ranges. Source-template ranges and the
/// independently intended target are intentionally absent from this selection.
pub(crate) fn select(
    snapshot: &Snapshot,
    bound: &Capture,
    report: &Report,
    index: usize,
    kind: Kind,
) -> Result<Selection, String> {
    validate(snapshot, bound, report)?;
    let junction = report
        .junctions
        .get(index)
        .ok_or("The selected actual junction is unavailable in this inspection.")?;
    // The native report checks each interval and its letters. Preserve contiguous
    // product-window order as well, so duplicated/reordered windows cannot move
    // the marked boundary to unrelated coordinates on repetitive DNA.
    let actual = snapshot.product_at(snapshot.plan.steps.len() - 1)?;
    let n = actual.sequence.len();
    let window = &junction.window;
    let first = window
        .parts
        .first()
        .ok_or("The junction window is empty.")?;
    let circular = actual.topology == collagene_core::Topology::Circular;
    let ordered = match window.parts.as_slice() {
        [_] => true,
        [left, right] => circular && left.end == n && right.start == 0,
        _ => false,
    };
    let boundary = first
        .start
        .checked_add(window.left_bases)
        .ok_or("Junction boundary coordinates overflow.")?;
    let boundary = if circular { boundary % n } else { boundary };
    if !ordered || boundary != junction.boundary {
        return Err(
            "The native junction window order no longer identifies its product boundary.".into(),
        );
    }
    let ranges = match kind {
        Kind::Window => &junction.window.parts,
        Kind::Homology => &junction.homology,
    };
    if ranges.is_empty() {
        return Err("This junction has no product homology interval to select.".into());
    }
    let n = snapshot
        .product_at(snapshot.plan.steps.len() - 1)?
        .sequence
        .len();
    let mut total = 0usize;
    let mut parts = Vec::with_capacity(ranges.len());
    for r in ranges {
        if r.start >= r.end || r.end > n {
            return Err("Junction selection lies outside the actual final product.".into());
        }
        total = total
            .checked_add(r.end - r.start)
            .filter(|length| *length <= n)
            .ok_or("Junction selection repeats more than one product turn.")?;
        if parts
            .iter()
            .any(|&(start, end)| r.start < end && start < r.end)
        {
            return Err("Junction selection contains repeated product bases.".into());
        }
        parts.push((r.start, r.end));
    }
    Ok(Selection {
        context_revision: bound.context_revision.clone(),
        product_revision: bound.product_revision.clone(),
        junction_index: index,
        kind,
        parts,
    })
}

/// Recheck immediately before opening a product canvas or applying highlights.
/// Neither matching DNA alone nor a previously selected numeric range is authority.
pub(crate) fn validate_selection(
    snapshot: &Snapshot,
    bound: &Capture,
    report: &Report,
    selected: &Selection,
) -> Result<(), String> {
    let fresh = select(
        snapshot,
        bound,
        report,
        selected.junction_index,
        selected.kind,
    )?;
    if selected != &fresh {
        return Err("The selected junction or its actual-product coordinates changed.".into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::{
        Molecule, Topology,
        chemistry::{EndTreatment, Terminus},
        verification::{DesignPlan, Operation, OrientedSource, PlanInput, Source, Status},
    };
    fn input(id: &str, text: &str) -> PlanInput {
        let molecule = Molecule::dna(id, text, Topology::Linear).unwrap();
        PlanInput {
            id: id.into(),
            revision: molecule.revision(),
            molecule,
        }
    }
    fn inspected(plan: DesignPlan, target: Molecule) -> (Snapshot, Capture, Report) {
        let options = Options::default();
        let report = assembly_junctions::inspect_plan(&plan, &target, &options).unwrap();
        let snapshot = Snapshot::new(plan, target, report.verification.clone());
        let capture = capture(&snapshot, options).unwrap();
        validate(&snapshot, &capture, &report).unwrap();
        (snapshot, capture, report)
    }
    fn overlap(circle: bool, wrong_target: bool) -> (Snapshot, Capture, Report) {
        let text = if wrong_target { "CCCCCC" } else { "AAACCC" };
        let target = Molecule::dna(
            "Independent target",
            text,
            if circle {
                Topology::Circular
            } else {
                Topology::Linear
            },
        )
        .unwrap();
        let (inputs, fragments, overlaps) = if circle {
            (
                vec![input("original", "AAACCCAAA")],
                vec![Source::Input("original".into())],
                vec![3],
            )
        } else {
            (
                vec![input("left", "AAAC"), input("right", "CCC")],
                vec![Source::Input("left".into()), Source::Input("right".into())],
                vec![1],
            )
        };
        let plan = DesignPlan {
            name: "Actual overlap".into(),
            target_revision: target.revision(),
            inputs,
            steps: vec![Operation::Overlap {
                fragments,
                overlaps,
                circular: circle,
            }],
        };
        inspected(plan, target)
    }
    fn ligation() -> (Snapshot, Capture, Report) {
        let target = Molecule::dna("Independent circle", "AAACCC", Topology::Circular).unwrap();
        let plan = DesignPlan {
            name: "Actual sealed circular ligation".into(),
            target_revision: target.revision(),
            inputs: vec![input("original", "AAACCC")],
            steps: vec![
                // An imported reference does not establish physical duplex
                // geometry or chemistry. Declare this original intact blunt
                // duplex and its 5′-OH/3′-OH ends as an explicit hypothesis.
                Operation::DeclareSourceTermini {
                    template: Source::Input("original".into()),
                    left: Terminus::PCR,
                    right: Terminus::PCR,
                },
                Operation::TreatEnds {
                    template: Source::Step(0),
                    treatment: EndTreatment::Phosphorylate,
                },
                Operation::Ligate {
                    fragments: vec![OrientedSource {
                        source: Source::Step(1),
                        reverse: false,
                    }],
                    circular: true,
                    allow_nicked: false,
                    options: Default::default(),
                },
            ],
        };
        inspected(plan, target)
    }
    #[test]
    fn actual_ligation_and_overlap_closure_keep_split_product_window_order() {
        for (snapshot, bound, report) in [ligation(), overlap(true, false)] {
            assert_eq!(
                report.verification.status,
                Status::ConsistentUnderAssumptions
            );
            let before = serde_json::to_vec(&snapshot.report).unwrap();
            let selected = select(&snapshot, &bound, &report, 0, Kind::Window).unwrap();
            assert_eq!(selected.parts, vec![(3, 6), (0, 3)]);
            assert_eq!(report.junctions[0].boundary, 0);
            assert_eq!(
                selected.product_revision,
                snapshot
                    .product_at(snapshot.plan.steps.len() - 1)
                    .unwrap()
                    .revision()
            );
            let actual = snapshot.product_at(snapshot.plan.steps.len() - 1).unwrap();
            let letters = selected
                .parts
                .iter()
                .map(|&(a, b)| &actual.sequence.as_str()[a..b])
                .collect::<String>();
            assert_eq!(letters, "CCCAAA");
            validate_selection(&snapshot, &bound, &report, &selected).unwrap();
            assert_eq!(serde_json::to_vec(&snapshot.report).unwrap(), before);
            if matches!(report.method, Some(assembly_junctions::Method::Overlap)) {
                assert_eq!(
                    select(&snapshot, &bound, &report, 0, Kind::Homology)
                        .unwrap()
                        .parts,
                    vec![(0, 3)]
                );
            } else {
                assert!(select(&snapshot, &bound, &report, 0, Kind::Homology).is_err());
            }
        }
    }
    #[test]
    fn wrong_independent_target_stays_invalid_while_actual_product_navigation_is_available() {
        let (snapshot, bound, report) = overlap(false, true);
        assert_eq!(report.verification.status, Status::Invalid);
        assert_eq!(snapshot.target.sequence.as_str(), "CCCCCC");
        assert_eq!(snapshot.product_at(0).unwrap().sequence.as_str(), "AAACCC");
        assert_eq!(
            select(&snapshot, &bound, &report, 0, Kind::Window)
                .unwrap()
                .parts,
            vec![(0, 6)]
        );
        assert_eq!(
            select(&snapshot, &bound, &report, 0, Kind::Homology)
                .unwrap()
                .parts,
            vec![(3, 4)]
        );
        assert_eq!(snapshot.report.status, Status::Invalid);
    }
    #[test]
    fn full_product_target_options_and_selected_coordinates_cannot_be_retargeted() {
        let (snapshot, bound, report) = overlap(true, false);
        let mut forged = bound.clone();
        forged.product_revision = "0".repeat(64);
        assert!(validate(&snapshot, &forged, &report).is_err());
        let mut forged = report.clone();
        forged.options.flank_bases = 29;
        assert!(validate(&snapshot, &bound, &forged).is_err());
        let mut forged = report.clone();
        forged.context_revision = "0".repeat(64);
        assert!(validate(&snapshot, &bound, &forged).is_err());
        let mut forged = report.clone();
        forged.junctions[0].window.parts[0].end = 1000;
        assert!(select(&snapshot, &bound, &forged, 0, Kind::Window).is_err());
        assert!(select(&snapshot, &bound, &report, usize::MAX, Kind::Window).is_err());
        let mut selected = select(&snapshot, &bound, &report, 0, Kind::Window).unwrap();
        selected.parts = vec![(0, 6)];
        assert!(validate_selection(&snapshot, &bound, &report, &selected).is_err());
        let mut selected = select(&snapshot, &bound, &report, 0, Kind::Window).unwrap();
        selected.context_revision = "0".repeat(64);
        assert!(validate_selection(&snapshot, &bound, &report, &selected).is_err());
        // Identical DNA letters do not authorize a changed full molecule.
        let mut product_changed = snapshot.report.clone();
        product_changed
            .products
            .last_mut()
            .unwrap()
            .metadata
            .insert("diagnostic change".into(), "changed".into());
        let changed = Snapshot::new(
            snapshot.plan.clone(),
            snapshot.target.clone(),
            product_changed,
        );
        assert!(validate(&changed, &bound, &report).is_err());
        let mut changed_target = snapshot.target.clone();
        changed_target.name = "Other independent target".into();
        let changed = Snapshot::new(
            snapshot.plan.clone(),
            changed_target,
            snapshot.report.clone(),
        );
        assert!(validate(&changed, &bound, &report).is_err());
    }
    #[test]
    fn current_product_bounds_and_dna_alphabet_are_checked_before_using_cached_receipts() {
        // A linear protein document is structurally valid, so this mutation
        // reaches the explicit DNA guard rather than the circular-protein ban.
        let (mut snapshot, _, _) = overlap(false, false);
        let product = snapshot.report.products.last_mut().unwrap();
        product.sequence = collagene_core::Sequence::new(
            collagene_core::Alphabet::Protein,
            product.sequence.as_str(),
        )
        .unwrap();
        // Even updating the cached full revision cannot turn the same symbolic
        // letters in a different alphabet into native junction DNA.
        snapshot.product_revision = Some(product.revision());
        let error = capture(&snapshot, Options::default()).unwrap_err();
        assert!(error.contains("DNA product"), "{error}");

        let (mut snapshot, _, _) = overlap(true, false);
        snapshot.report.products.last_mut().unwrap().name = "x".repeat(1025);
        let error = capture(&snapshot, Options::default()).unwrap_err();
        assert!(error.contains("supported size limit"), "{error}");
    }

    #[test]
    fn a_leading_preparation_never_stands_in_for_missing_final_product() {
        let (snapshot, _, _) = ligation();
        let mut report = snapshot.report.clone();
        report.products.pop();
        report.status = Status::Incomplete;
        let expected_prefix = snapshot.plan.steps.len() - 1;
        let prefix = Snapshot::new(snapshot.plan.clone(), snapshot.target.clone(), report);
        assert_eq!(prefix.prefix_available, expected_prefix);
        assert!(prefix.product_revision.is_none());
        assert!(capture(&prefix, Options::default()).is_err());
    }
}
