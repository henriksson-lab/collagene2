#[path = "../src/history/history_view/graph.rs"]
mod graph;
use collagene_core::documents::history::HistorySummary;
use graph::*;
fn event(id: &str, revision: &str, parents: &[&str]) -> HistorySummary {
    HistorySummary {
        event_id: id.into(),
        revision: revision.into(),
        operation: id.into(),
        parents: parents.iter().map(|s| s.to_string()).collect(),
        current: false,
        sequence_length: 8,
    }
}
#[test]
fn branch_and_merge_draw_true_event_edges_without_collapsing_content() {
    let events = vec![
        event("root", "A", &[]),
        event("left", "B", &["root"]),
        event("right", "A", &["root"]),
        event("merge", "C", &["left", "right"]),
    ];
    let index = GraphIndex::new(&events, "merge").unwrap();
    let page = index.page(0, &[], &[]);
    assert_eq!(
        page.nodes
            .iter()
            .map(|n| n.event_id.as_str())
            .collect::<Vec<_>>(),
        ["merge", "right", "left", "root"]
    );
    let node = |id: &str| page.nodes.iter().find(|n| n.event_id == id).unwrap();
    assert_eq!(node("root").same_content_events, 2);
    assert_eq!(node("right").same_content_events, 2);
    assert!(node("root").branch_point);
    assert!(!node("merge").branch_point);
    assert_eq!(node("left").lane, node("root").lane);
    assert_ne!(node("right").lane, node("left").lane);
    let mut links = page
        .edges
        .iter()
        .map(|e| (e.child.as_str(), e.parent.as_str()))
        .collect::<Vec<_>>();
    links.sort();
    assert_eq!(
        links,
        [
            ("left", "root"),
            ("merge", "left"),
            ("merge", "right"),
            ("right", "root")
        ]
    );
    let mut geometry = page
        .edges
        .iter()
        .map(|e| {
            (
                e.child.as_str(),
                e.parent.as_str(),
                e.from_row,
                e.from_lane,
                e.to,
            )
        })
        .collect::<Vec<_>>();
    geometry.sort();
    assert_eq!(
        geometry,
        [
            ("left", "root", 2, 0, Some((3, 0))),
            ("merge", "left", 0, 0, Some((2, 0))),
            ("merge", "right", 0, 0, Some((1, 1))),
            ("right", "root", 1, 1, Some((3, 0)))
        ]
    );
    assert!(page.edges.iter().all(|e| e.to.is_some()));
    assert!(page.nodes.iter().all(|n| n.ancestor));
}
#[test]
fn navigation_is_not_inferred_from_ancestry_and_duplicate_visits_keep_next_distance() {
    let events = vec![
        event("root", "A", &[]),
        event("left", "B", &["root"]),
        event("right", "C", &["root"]),
    ];
    let index = GraphIndex::new(&events, "left").unwrap();
    let page = index.page(
        0,
        &["root".into(), "right".into(), "root".into()],
        &["right".into()],
    );
    let right = page.nodes.iter().find(|n| n.event_id == "right").unwrap();
    assert!(!right.ancestor);
    assert_eq!(right.undo_distance, Some(2));
    assert_eq!(right.redo_distance, Some(1));
    let root = page.nodes.iter().find(|n| n.event_id == "root").unwrap();
    assert!(root.ancestor);
    assert_eq!(root.undo_distance, Some(1));
    assert!(
        page.nodes
            .iter()
            .find(|n| n.event_id == "left")
            .unwrap()
            .current
    );
}
#[test]
fn unordered_serialization_and_external_parents_have_explicit_display_meaning() {
    let events = vec![
        event("child", "C", &["parent", "other-document-event"]),
        event("root", "A", &[]),
        event("parent", "B", &["root"]),
    ];
    let index = GraphIndex::new(&events, "child").unwrap();
    let page = index.page(0, &[], &[]);
    assert_eq!(
        page.nodes
            .iter()
            .map(|n| n.event_id.as_str())
            .collect::<Vec<_>>(),
        ["child", "parent", "root"]
    );
    let external = page
        .edges
        .iter()
        .find(|e| e.parent == "other-document-event")
        .unwrap();
    assert!(!external.parent_in_document);
    assert!(external.to.is_none());
    assert_eq!(index.event_at(0), Some("child"));
    assert_eq!(index.display_position("root"), Some(2));
}
#[test]
fn pages_disclose_off_page_parents_and_reach_every_event() {
    let events = (0..105)
        .map(|i| {
            event(
                &format!("e{i}"),
                &format!("r{i}"),
                &if i == 0 {
                    vec![]
                } else {
                    vec![format!("e{}", i - 1)]
                }
                .iter()
                .map(String::as_str)
                .collect::<Vec<_>>(),
            )
        })
        .collect::<Vec<_>>();
    let index = GraphIndex::new(&events, "e104").unwrap();
    let first = index.page(0, &[], &[]);
    let second = index.page(1, &[], &[]);
    assert_eq!(first.nodes.len(), 100);
    assert_eq!(second.nodes.len(), 5);
    assert_eq!(
        first
            .edges
            .iter()
            .filter(|e| e.to.is_none() && e.parent_in_document)
            .count(),
        1
    );
    assert_eq!(first.nodes.last().unwrap().event_id, "e5");
    assert_eq!(second.nodes.first().unwrap().event_id, "e4");
    let all = first
        .nodes
        .iter()
        .chain(&second.nodes)
        .map(|n| n.event_id.clone())
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(all.len(), 105);
    assert!(index.page(2, &[], &[]).nodes.is_empty());
}
#[test]
fn graph_rejects_cycles_duplicate_identity_and_over_budget_metadata() {
    assert!(GraphIndex::new(&[event("a", "A", &["b"]), event("b", "B", &["a"])], "a").is_err());
    assert!(GraphIndex::new(&[event("a", "A", &[]), event("a", "B", &[])], "a").is_err());
    assert!(GraphIndex::new(&[event("a", "A", &["outside", "outside"])], "a").is_err());
    assert!(GraphIndex::new(&[event("a", "A", &[])], "missing").is_err());
    let too_many = vec![event("a", "A", &[]); MAX_GRAPH_EVENTS + 1];
    assert!(
        GraphIndex::new(&too_many, "a")
            .unwrap_err()
            .contains("20,000")
    );
}
#[test]
fn edge_drawing_limit_does_not_drop_event_nodes_or_claim_all_edges_are_drawn() {
    // 100 visible roots each reference distinct external parents; every event remains.
    let events = (0..100)
        .map(|i| {
            event(
                &format!("e{i}"),
                "same",
                &(0..30)
                    .map(|p| format!("external{p}"))
                    .collect::<Vec<_>>()
                    .iter()
                    .map(String::as_str)
                    .collect::<Vec<_>>(),
            )
        })
        .collect::<Vec<_>>();
    let page = GraphIndex::new(&events, "e99").unwrap().page(0, &[], &[]);
    assert_eq!(page.nodes.len(), 100);
    assert_eq!(page.edges.len(), MAX_PAGE_EDGES);
    assert_eq!(page.omitted_edges, 1000);
    assert_eq!(page.lane_count, 100);
    assert_eq!(page.width(), 2024);
    assert_eq!(page.height(), 4624);
}
