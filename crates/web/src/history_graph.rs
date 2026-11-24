//! Bounded event-DAG layout. Event identity is never collapsed into sequence content.
use collagene_core::history::HistorySummary;
use std::collections::{BTreeMap, BTreeSet};

pub const PAGE_SIZE: usize = 100;
pub const MAX_GRAPH_EVENTS: usize = 20_000;
pub const MAX_GRAPH_PARENTS: usize = 100_000;
pub const MAX_PAGE_EDGES: usize = 2_000;
pub const ROW_HEIGHT: usize = 46;

#[derive(Debug, Clone)]
pub struct GraphNode {
    pub event_index: usize,
    pub event_id: String,
    pub row: usize,
    pub lane: usize,
    pub current: bool,
    pub ancestor: bool,
    pub branch_point: bool,
    pub same_content_events: usize,
    /// Distance in the actual navigation stack, not inferred ancestry. One means next.
    pub undo_distance: Option<usize>,
    pub redo_distance: Option<usize>,
}
#[derive(Debug, Clone)]
pub struct GraphEdge {
    pub child: String,
    pub parent: String,
    pub from_row: usize,
    pub from_lane: usize,
    pub to: Option<(usize, usize)>,
    /// A local event outside this page is distinct from an external document parent.
    pub parent_in_document: bool,
}
#[derive(Debug, Clone)]
pub struct GraphPage {
    pub nodes: Vec<GraphNode>,
    pub edges: Vec<GraphEdge>,
    pub omitted_edges: usize,
    pub lane_count: usize,
}
impl GraphPage {
    pub fn width(&self) -> usize {
        24 + self.lane_count.max(1) * 20
    }
    pub fn height(&self) -> usize {
        24 + self.nodes.len().max(1) * ROW_HEIGHT
    }
}
#[derive(Debug, Clone)]
pub struct GraphIndex {
    order: Vec<usize>,
    ids: BTreeMap<String, usize>,
    display_positions: BTreeMap<String, usize>,
    lanes: Vec<usize>,
    parents: Vec<Vec<String>>,
    ancestor: Vec<bool>,
    current: Vec<bool>,
    branch_point: Vec<bool>,
    content_count: Vec<usize>,
    event_ids: Vec<String>,
}
impl GraphIndex {
    pub fn new(events: &[HistorySummary], head: &str) -> Result<Self, String> {
        if events.len() > MAX_GRAPH_EVENTS {
            return Err("Graph exceeds the 20,000-event display index bound; the paginated event list remains available.".into());
        }
        let mut total = 0usize;
        for event in events {
            total=total.checked_add(event.parents.len()).filter(|n|*n<=MAX_GRAPH_PARENTS).ok_or("Graph exceeds the 100,000-parent display index bound; the paginated event list remains available.")?;
        }
        let mut ids = BTreeMap::new();
        for (i, event) in events.iter().enumerate() {
            if ids.insert(event.event_id.clone(), i).is_some() {
                return Err(
                    "Duplicate event identities cannot be drawn as separate history nodes.".into(),
                );
            }
        }
        if !events.is_empty() && !ids.contains_key(head) {
            return Err("Current event is absent from this document's retained history.".into());
        }
        let mut indegrees = vec![0usize; events.len()];
        let mut children = vec![Vec::new(); events.len()];
        for (i, event) in events.iter().enumerate() {
            let mut seen = BTreeSet::new();
            for parent in &event.parents {
                if !seen.insert(parent) {
                    return Err(
                        "Duplicate parent references cannot be drawn as distinct history edges."
                            .into(),
                    );
                }
                if let Some(&p) = ids.get(parent) {
                    indegrees[i] += 1;
                    children[p].push(i);
                }
            }
        }
        // Imported archives need not serialize their event records in topological order.
        let mut ready = (0..events.len())
            .filter(|&i| indegrees[i] == 0)
            .collect::<BTreeSet<_>>();
        let mut chronological = Vec::with_capacity(events.len());
        while let Some(i) = ready.pop_first() {
            chronological.push(i);
            for &child in &children[i] {
                indegrees[child] -= 1;
                if indegrees[child] == 0 {
                    ready.insert(child);
                }
            }
        }
        if chronological.len() != events.len() {
            return Err("Cyclic event parents cannot be drawn as a revision DAG.".into());
        }
        let mut lanes = vec![0; events.len()];
        let mut next_lane = 0;
        let mut assigned_children = vec![0usize; events.len()];
        for &i in &chronological {
            let local = events[i]
                .parents
                .iter()
                .filter_map(|p| ids.get(p).copied())
                .collect::<Vec<_>>();
            lanes[i] = if local.is_empty() {
                let at = next_lane;
                next_lane += 1;
                at
            } else if local.len() > 1 {
                local.iter().map(|&p| lanes[p]).min().unwrap_or(0)
            } else {
                let parent = local[0];
                if assigned_children[parent] == 0 {
                    lanes[parent]
                } else {
                    let at = next_lane;
                    next_lane += 1;
                    at
                }
            };
            for p in local {
                assigned_children[p] += 1;
            }
        }
        let mut ancestor = vec![false; events.len()];
        let mut pending = ids.get(head).copied().into_iter().collect::<Vec<_>>();
        while let Some(i) = pending.pop() {
            if ancestor[i] {
                continue;
            }
            ancestor[i] = true;
            pending.extend(events[i].parents.iter().filter_map(|p| ids.get(p).copied()));
        }
        let mut counts = BTreeMap::<&str, usize>::new();
        for event in events {
            *counts.entry(&event.revision).or_default() += 1;
        }
        let order = chronological.into_iter().rev().collect::<Vec<_>>();
        let display_positions = order
            .iter()
            .enumerate()
            .map(|(position, &i)| (events[i].event_id.clone(), position))
            .collect();
        Ok(Self {
            order,
            ids,
            display_positions,
            lanes,
            parents: events.iter().map(|e| e.parents.clone()).collect(),
            ancestor,
            current: events.iter().map(|e| e.event_id == head).collect(),
            branch_point: children.iter().map(|c| c.len() > 1).collect(),
            content_count: events.iter().map(|e| counts[e.revision.as_str()]).collect(),
            event_ids: events.iter().map(|e| e.event_id.clone()).collect(),
        })
    }
    pub fn display_position(&self, id: &str) -> Option<usize> {
        self.display_positions.get(id).copied()
    }
    pub fn event_at(&self, position: usize) -> Option<&str> {
        self.order
            .get(position)
            .map(|&i| self.event_ids[i].as_str())
    }
    pub fn display_indices(&self, page: usize) -> Vec<usize> {
        self.order
            .iter()
            .skip(page.saturating_mul(PAGE_SIZE))
            .take(PAGE_SIZE)
            .copied()
            .collect()
    }
    pub fn page(&self, page: usize, undo: &[String], redo: &[String]) -> GraphPage {
        let indices = self.display_indices(page);
        let lanes = indices
            .iter()
            .map(|&i| self.lanes[i])
            .collect::<BTreeSet<_>>();
        let compact = lanes
            .iter()
            .enumerate()
            .map(|(at, &lane)| (lane, at))
            .collect::<BTreeMap<_, _>>();
        let undo = undo
            .iter()
            .rev()
            .enumerate()
            .fold(BTreeMap::new(), |mut out, (i, id)| {
                out.entry(id.as_str()).or_insert(i + 1);
                out
            });
        let redo = redo
            .iter()
            .rev()
            .enumerate()
            .fold(BTreeMap::new(), |mut out, (i, id)| {
                out.entry(id.as_str()).or_insert(i + 1);
                out
            });
        let nodes = indices
            .iter()
            .enumerate()
            .map(|(row, &i)| GraphNode {
                event_index: i,
                event_id: self.event_ids[i].clone(),
                row,
                lane: compact[&self.lanes[i]],
                current: self.current[i],
                ancestor: self.ancestor[i],
                branch_point: self.branch_point[i],
                same_content_events: self.content_count[i],
                undo_distance: undo.get(self.event_ids[i].as_str()).copied(),
                redo_distance: redo.get(self.event_ids[i].as_str()).copied(),
            })
            .collect::<Vec<_>>();
        let visible = nodes
            .iter()
            .map(|n| (n.event_id.as_str(), (n.row, n.lane)))
            .collect::<BTreeMap<_, _>>();
        let mut edges = Vec::new();
        let mut omitted_edges = 0;
        for node in &nodes {
            for parent in &self.parents[node.event_index] {
                if edges.len() == MAX_PAGE_EDGES {
                    omitted_edges += 1;
                    continue;
                }
                edges.push(GraphEdge {
                    child: node.event_id.clone(),
                    parent: parent.clone(),
                    from_row: node.row,
                    from_lane: node.lane,
                    to: visible.get(parent.as_str()).copied(),
                    parent_in_document: self.ids.contains_key(parent),
                });
            }
        }
        GraphPage {
            nodes,
            edges,
            omitted_edges,
            lane_count: lanes.len(),
        }
    }
}
