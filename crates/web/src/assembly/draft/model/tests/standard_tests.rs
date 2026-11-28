use super::*;
use crate::assembly::draft::assembly_standard_model::StandardDraft;
use collagene_assembly_standards::assessment::Outcome;
use collagene_core::{
    chemistry::Terminus,
    cloning::{
        assembly::{
            recipe::CleavagePolicy,
            standards::{
                self,
                declaration::{DeclaredReviewState, ProfileRef},
            },
        },
        ligation::recipe,
    },
    documents::workspace::Workspace,
};

fn fixture() -> (Workspace, Draft, Molecule) {
    let frame = Workspace::from_documents(
        "source format fixture",
        vec![
            Molecule::dna("Promoter reference", "ACACAC", Topology::Linear).unwrap(),
            Molecule::dna("UTR reference", "CACACA", Topology::Linear).unwrap(),
        ],
        0,
    )
    .unwrap();
    let mut draft = Draft::default();
    draft.circular = false;
    draft.methods.mode = Mode::Ligation;
    for index in 0..2 {
        let document = &frame.documents[index];
        let molecule = frame.revisions[&document.head].clone();
        draft.slots[index].captured = Some(Arc::new(CapturedSource {
            source: SourceInfo {
                binding: collagene_core::material::lineage::EventBinding {
                    document_id: document.id.clone(),
                    event_id: document.head_event.clone().unwrap(),
                    revision: document.head.clone(),
                },
                name: molecule.name.clone(),
                length: 6,
                alphabet: Alphabet::Dna,
                topology: Topology::Linear,
            },
            molecule,
        }));
        draft.slots[index].selector = assembly_intent::Selector::ReferenceRange {
            start: 0,
            end: 6,
            reverse: false,
        };
        draft.methods.rows[index].preparation = OverlapPreparation::DeclaredDuplex;
        draft.methods.rows[index].left_terminus = Terminus::CLEAVED;
        draft.methods.rows[index].right_terminus = Terminus::CLEAVED;
    }
    draft.sync_rows();
    draft.methods.source_roles = Some(draft.methods.role_ids(draft.slots.len()).unwrap());
    let profile = standards::catalogue()
        .unwrap()
        .into_iter()
        .find(|p| p.profile_id == collagene_assembly_standards::profile::RFC1000_L0_ID)
        .unwrap();
    let mut standard = StandardDraft::capture(
        &draft,
        ProfileRef {
            profile_id: profile.profile_id,
            profile_identity: profile.profile_identity,
        },
    )
    .unwrap();
    for (row, slot) in standard.rows.iter_mut().zip(["promoter", "five-prime-utr"]) {
        row.slot = slot.into();
        row.body_start = "0".into();
        row.body_end = "6".into();
    }
    draft.methods.standard = Some(standard);
    draft.validate().unwrap();
    let target = Molecule::dna(
        "Independent literal target",
        "ACACACCACACA",
        Topology::Linear,
    )
    .unwrap();
    (frame, draft, target)
}
fn preview(draft: &Draft, target: &Molecule) -> recipe::Preview {
    let inputs = draft
        .slots
        .iter()
        .map(|slot| slot.captured.as_ref().unwrap().molecule.clone())
        .collect();
    let (inputs, request) = crate::cloning::preparation::cloning_preparations::ligation_request(
        &draft.methods,
        inputs,
        target,
    )
    .unwrap();
    recipe::preview(&inputs, target, &request, CleavagePolicy::AssumedComplete).unwrap()
}
#[test]
fn exact_standard_is_in_plan_before_hash_and_advisory_failure_preserves_actual_physical_products() {
    let (_, draft, target) = fixture();
    let good = preview(&draft, &target);
    let declaration = draft.methods.standard_declaration().unwrap();
    assert_eq!(good.compiled.plan.standard, declaration);
    assert_eq!(
        good.verification.plan_revision,
        good.compiled.plan.revision()
    );
    let DeclaredReviewState::Assessed(report) =
        &good.verification.standard_review.as_ref().unwrap().state
    else {
        panic!("fresh native advice required");
    };
    assert_eq!(report.assessment.outcome, Outcome::Conforms);
    let mut bad = draft.clone();
    bad.methods.standard.as_mut().unwrap().rows[0].slot = "terminator".into();
    let failed = preview(&bad, &target);
    let DeclaredReviewState::Assessed(report) =
        &failed.verification.standard_review.as_ref().unwrap().state
    else {
        panic!("native violation evidence required");
    };
    assert_eq!(report.assessment.outcome, Outcome::Violations);
    assert_eq!(failed.verification.status, good.verification.status);
    let a = good.verification.final_product().unwrap();
    let b = failed.verification.final_product().unwrap();
    assert_eq!(
        (&a.sequence, a.topology, &a.features),
        (&b.sequence, b.topology, &b.features)
    );
    assert_eq!(a.sequence, target.sequence);
    assert_eq!(
        serde_json::to_value(&good.product_physical).unwrap(),
        serde_json::to_value(&failed.product_physical).unwrap()
    );
    assert_ne!(
        good.compiled.plan.revision(),
        failed.compiled.plan.revision()
    );
}
#[test]
fn role_reorder_changes_declared_grammar_order_without_guessing_new_body_geometry() {
    let (_, mut draft, _) = fixture();
    let source_bindings = draft.methods.standard.clone().unwrap();
    draft.slots.swap(0, 1);
    draft.methods.swap_rows(0, 1);
    draft.sync_rows();
    source_bindings.validate_sources(&draft).unwrap();
    let declaration = draft.methods.standard_declaration().unwrap().unwrap();
    assert_eq!(declaration.parts[0].slot, "five-prime-utr");
    assert_eq!(declaration.parts[1].slot, "promoter");
    assert!(
        declaration
            .parts
            .iter()
            .all(|part| part.body_start == 0 && part.body_end == 6)
    );
}
#[test]
fn stale_exact_event_selector_and_preparation_require_explicit_rebinding_but_remain_editable() {
    let (_, draft, _) = fixture();
    let bound = draft.methods.standard.clone().unwrap();
    let mut changed = draft.clone();
    let mut history = collagene_core::documents::history::DocumentHistory::new(
        "document-1",
        changed.slots[0].captured.as_ref().unwrap().molecule.clone(),
    )
    .unwrap();
    use collagene_core::documents::reference_edit::Operation;
    history
        .commit(
            vec![Operation::Rename {
                name: "temporary".into(),
            }],
            "rename",
        )
        .unwrap();
    history
        .commit(
            vec![Operation::Rename {
                name: "Promoter reference".into(),
            }],
            "rename back",
        )
        .unwrap();
    let mut capture = changed.slots[0].captured.as_ref().unwrap().as_ref().clone();
    assert_eq!(capture.source.binding.revision, history.molecule.revision());
    assert_ne!(capture.source.binding.event_id, history.binding().event_id);
    capture.source.binding = history.binding();
    changed.slots[0].captured = Some(Arc::new(capture));
    assert!(bound.validate_sources(&changed).is_err());
    assert!(
        changed.validate().is_ok(),
        "stale editable declarations are retained, not silently cleared"
    );
    let mut changed = draft.clone();
    changed.slots[0].selector = assembly_intent::Selector::ReferenceRange {
        start: 1,
        end: 6,
        reverse: false,
    };
    assert!(bound.validate_sources(&changed).is_err());
    let mut changed = draft.clone();
    changed.methods.rows[0].preparation = OverlapPreparation::WholeLinear;
    assert!(bound.validate_sources(&changed).is_err());
    let fresh = StandardDraft::capture(&draft, bound.profile.clone()).unwrap();
    assert!(
        fresh
            .rows
            .iter()
            .all(|r| r.slot.is_empty() && r.body_start.is_empty() && r.body_end.is_empty())
    );
    assert!(
        fresh.declaration(&draft.methods).is_err(),
        "body bounds are explicit, never guessed from a profile or feature"
    );
}
#[test]
fn current_named_design_roundtrip_retains_exact_events_and_profile_but_no_review_authority() {
    let (frame, draft, _) = fixture();
    let bytes = encode(&draft).unwrap();
    let restored = decode(&bytes).unwrap();
    assert_eq!(restored.methods.standard, draft.methods.standard);
    let (saved, entry) = crate::assembly::library::assembly_saved_designs_model::save(
        &frame,
        "Standard source assembly",
        &draft,
    )
    .unwrap();
    let reopened =
        crate::assembly::library::assembly_saved_designs_model::load(&saved, &entry.id).unwrap();
    assert_eq!(reopened.methods.standard, draft.methods.standard);
    let raw = serde_json::to_value(&reopened).unwrap();
    assert!(raw["methods"]["standard"].is_object());
    assert!(raw.get("standard_review").is_none());
}
#[test]
fn html_and_pdf_bind_fresh_advisory_scope_and_same_canonical_plan() {
    let (_, draft, target) = fixture();
    let result = preview(&draft, &target);
    let plan = &result.compiled.plan;
    let html = collagene_core::verification::output::export_protocol_html(plan, &target);
    let review = result.verification.standard_review.as_ref().unwrap();
    assert!(html.contains(&review.declaration_revision));
    assert!(html.contains(&review.review_revision));
    assert!(
        html.contains(
            &draft
                .methods
                .standard
                .as_ref()
                .unwrap()
                .profile
                .profile_identity
        )
    );
    assert!(html.contains("advisory format review"));
    assert!(!html.contains("full standard verified"));
    let pdf = collagene_core::output::pdf::export_protocol_pdf(plan, &target).unwrap();
    assert!(pdf.bytes.starts_with(b"%PDF-"));
    assert_eq!(pdf.plan_revision, plan.revision());
    assert_eq!(pdf.target_revision, target.revision());
    assert_eq!(pdf.status, result.verification.status);
}
