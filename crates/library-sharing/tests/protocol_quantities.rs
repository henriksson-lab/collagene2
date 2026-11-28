#[path = "../../core/tests/common/protocol_quantities.rs"]
mod fixture;
use collagene_core::{
    documents::workspace, material::lineage as run_lineage, output::quantities as q,
};
use collagene_library_sharing::{self as sharing, Selection};
#[test]
fn captured_source_bound_quantities_are_omitted_from_ordinary_plasmid_sharing_with_counts() {
    let (w, r) = fixture::fixture();
    let captured = q::capture(&w, &r).unwrap();
    let w = captured.workspace;
    let before = workspace::save_workspace(&w).unwrap();
    let run_id = match &r.plan_source {
        q::PlanSource::SavedRun { run_id } => run_id,
        _ => unreachable!(),
    };
    let run = run_lineage::load_run(&w, run_id).unwrap();
    let p = run.products.last().unwrap();
    let selection = Selection {
        document_id: p.document_id.clone(),
        revision: p.revision.clone(),
        head_event: Some(p.event_id.clone()),
    };
    // Rendering a fresh combined protocol must not attach a scientific status
    // or disclose the documentary quantities through ordinary plasmid sharing.
    for format in [q::ExportFormat::ProtocolHtml, q::ExportFormat::ProtocolPdf] {
        let rendered = q::export(&w, &captured.attachment_hash, format).unwrap();
        assert!(rendered.molecular.is_some());
    }
    let projected = sharing::document_closure(&w, &selection).unwrap();
    assert!(
        !projected
            .inventory
            .attachments
            .contains(&captured.attachment_hash)
    );
    assert!(projected.inventory.omitted.attachments >= 1);
    assert!(
        projected
            .workspace
            .attachments
            .values()
            .all(|a| a.kind != q::KIND)
    );
    for index in [2, 3, 4] {
        assert!(
            projected
                .workspace
                .documents
                .iter()
                .all(|d| d.id != w.documents[index].id)
        );
    }
    assert!(
        !serde_json::to_string(&projected.workspace)
            .unwrap()
            .contains("PRIVATE EXPLICIT QUANTITY PLAN")
    );
    let restored = workspace::load_workspace(&projected.archive).unwrap();
    assert!(
        run_lineage::replay_run(&restored, run_id)
            .unwrap()
            .products_match
    );
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
}
#[test]
fn unknown_quantity_extension_kind_is_not_silently_omitted_as_a_generic_prefix() {
    let (w, r) = fixture::fixture();
    let c = q::capture(&w, &r).unwrap();
    let mut w = c.workspace;
    w.attachments.get_mut(&c.attachment_hash).unwrap().kind = "protocol-quantities-v2".into();
    let d = &w.documents[0];
    let selected = Selection {
        document_id: d.id.clone(),
        revision: d.head.clone(),
        head_event: d.head_event.clone(),
    };
    assert!(matches!(
        sharing::document_closure(&w, &selected),
        Err(sharing::Error::UnsupportedAttachment(_))
    ));
}
