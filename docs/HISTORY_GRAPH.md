# Retained event history graph

The **History** panel displays the recorded document event DAG beside an accessible graph navigator and the existing paginated revision cards. It does not infer causal edges from sequence similarity, sequence equality, filenames, candidate labels, or the order in which documents were opened. Every drawn arrow points from an event to one of its exact recorded parent event IDs.

The graph keeps event identity separate from a molecule content revision. Returning to a previous molecule state can therefore produce two distinct event nodes with the same revision hash. A badge identifies shared molecule revisions without collapsing those events. The current event, its recorded local ancestry, other retained branches/roots, and branch points have separate labels. Parent-first topological ordering is derived from the recorded local edges; serialization order breaks ties. The graph does not claim timestamps or experimental chronology.

Undo and Redo badges come from the controller's exact event navigation stacks, with distance one identifying the next action. These stacks are distinct from biological ancestry: restoring a retained branch adds the prior head to navigation history but does not rewrite the branch's parent edges. The read-only `DocumentHistory::navigation_event_ids()` helper returns snapshots of these stacks, oldest entry first and next operation last.

## Navigation and forks

Click a graph circle or its labeled button to inspect the event without modifying the molecule. Up/Down, Home and End move graph inspection; selection crosses page boundaries and retains keyboard focus. **Show current event** returns to the active head. Paste a full retained event ID into **Find event** to inspect it without modifying the active molecule; an external parent can be found after opening its source document. **Show parent event** navigates a parent retained in the current document. The existing Restore/Open actions and the selected-event actions invoke the same application transactions and capacity checks.

`DocumentHistory::fork_event(source_event, new_document_id)` opens an immutable source event as a distinct document whose first event records the exact source-event parent. The source event must remain retained in the workspace; normal workspace validation rejects a dangling parent. Forking copies no verification authority and creates no Undo/Redo entries for the new document. The workspace archive retains this cross-document edge and deduplicates shared molecule revisions. The new document can then be edited normally, preserving its provenance chain.

Dashed graph continuations distinguish parents elsewhere on the current graph's pages from parents outside the current document. External parent display metadata comes from other retained workspace documents. An unavailable display name is explicitly labeled unavailable, never replaced by an inferred source. Full parent IDs remain inspectable. External labels are bounded display text; event IDs and revision IDs retain their full values.

Saved candidate targets can match a node's molecule revision. This is labeled a content match, not proof that the event created or verified that candidate. Candidate snapshots remain proposals for separate review and deterministic verification against exact active inputs and an independently selected target. History operations do not execute recorded operation labels or establish a scientific verification result.

## Display bounds

The graph index supports up to 20,000 events and 100,000 recorded parent references, consistent with the workspace-level event bound and a separate visualization work bound. Exceeding either display-index bound leaves the paginated revision cards available. Invalid duplicate identities, duplicate parent references or local cycles cannot be rendered as a valid DAG.

Each page contains at most 100 event nodes and 2,000 drawn parent edges. Additional parent edges are counted and explicitly disclosed; all parent IDs remain in the event details. Visible logical lanes are compacted per page, so even deep retained histories do not require a chromosome-size SVG. Branch-heavy pages can scroll horizontally inside their bounded graph viewport. Every retained event remains reachable through pagination; rendering is not a first-100-event truncation.

The graph is a Rust/Leptos SVG and DOM view. It adds no JavaScript dependency, external scientific runtime or new package dependency. Native layout tests use independently specified branch/merge coordinates, duplicate-content events, navigation-versus-ancestry cases, unordered records, external/page continuations and edge-budget disclosure. Native history tests cover exact fork/navigation persistence. `scripts/history_graph_browser_qa.py` checks the actual web actions, equal-content identities, retained forks after reload, keyboard focus and page transitions, and mobile overflow bounds.
