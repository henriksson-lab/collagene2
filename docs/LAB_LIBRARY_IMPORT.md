# Importing a fetched plasmid closure into a local library

`collagene_library_sharing::import::merge_document_closure` proposes an atomic
native graph merge. It preserves the existing local library, unpublished branches
and exact scientific ancestry; it does not replace the local workspace with a
fetched archive. Inputs are borrowed and remain unchanged on success and error.
The caller receives a new validated workspace only after all checks complete.

```rust,ignore
let result = merge_document_closure(
    &local_workspace, &loaded_fetched_closure, &selected_native_pin,
    &ImportOptions::default(),
)?;
```

The selected pin is a native document ID, molecule revision and optional native
head event. It is separate from any database document/revision ID or remote URL.
A byte-valid archive and an authenticated download do not confer scientific
material authority; ordinary saved-origin reconstruction still applies.

## Exact graph union

Both workspaces first pass the existing native workspace validator. The fetched
workspace must be the exact supported selected-document publication closure:
independently derived required document/revision/event/candidate/attachment sets
must match, every projected dependency head/event must match, selected activity
must match, and fetched navigation stacks must be empty. A full remote workspace
or unrelated private sidecar cannot be imported by pretending it is one plasmid.
Closure construction remains the separate publication helper's responsibility.

The merge treats revisions, events, candidates and attachments as immutable:

- Same identity and exact content is deduplicated.
- Same identity with different content returns a typed `Collision` identifying
  its kind and identity. It never chooses the newest version or rewrites hashes.
- Molecules, history parents, event operation records, candidate declarations and
  attachment bytes/metadata retain their native identities and contents.
- Missing native document records are appended without renaming their IDs.
- Existing **dependency** document records retain their complete local head/event
  and undo/redo navigation, even when a fetched run requires an older source event.
  That exact older event/revision is retained by the immutable union and remains
  available for fresh run reconstruction.

Local document order, workspace name, active candidate and unrelated local
attachments/proposals remain. By default the imported selected native document
becomes active; `activate_selected: false` keeps the local active document.
Projected counts, aggregate bytes and metadata are bounded before large model
clones, and the complete proposed result passes native workspace validation.
No scalar `.cg2` fallback, event relabeling, or source-target substitution occurs.

## Divergent selected heads

The default `SelectedPolicy::Strict` accepts a missing selected document or an
exactly matching local head/event. A different selected head **or event** yields
`SelectedConflict { local, incoming }`. Sequence equality does not resolve this:
two events can describe different branches even when the complete molecule is
identical. The caller keeps both local and fetched work available for a deliberate
choice; no remote state is automatically preferred.

Two explicit policies are available:

- `FastForward { expected_local }` requires the exact current local document/head/
  event and proves its event is an ancestor of the fetched event through retained
  parent edges. Eventless records and same-sequence unrelated events do not prove
  ancestry. A changed local pin yields `StaleLocal`; a non-descendant yields
  `NotFastForward`.
- `Checkout { expected_local }` allows an explicit branch switch after the exact
  expected local pin is checked. It does not assert an ancestry relationship.

Both transitions clear the selected document's active undo/redo stacks and return
`displaced_selected: Some(HeadBookmark)`. The bookmark contains the entire prior
local document navigation record, including its exact previous head/event and
undo/redo stacks, plus the imported destination pin. All old revisions/events
also remain in the resulting native graph. Dependency document navigation is
never cleared or switched by these policies.

**The caller must durably retain this bookmark before committing the returned
workspace.** A checkout result is a proposal, not a self-persisting transaction.
Persist the bookmark with the local library's branch/bookkeeping record so the
old branch can be recovered after restart; do not merely display it or keep it
in a transient signal. If persistence fails, do not commit the workspace. The
bookmark is ordinary branch navigation, not a saved material receipt.

If the caller performs the work asynchronously or keeps a preview, it must bind
the result to the entire local workspace snapshot and recheck it at commit time.
`expected_local` protects the selected head; it is not permission to overwrite
unrelated local changes made after an earlier snapshot. Build from the current
workspace and commit atomically with its existing persistence/conflict guard.

## Offline work and remotes

This module is pure Rust/native-workspace logic and is usable offline, including
in the WASM client and a future Tauri local library. It performs no HTTP request,
database mutation, filesystem write or scientific replay. Local plasmid analysis
and local folder/branch work require no remote database authority.

Remote choices, database base revisions, sync cursors and local branch bookmarks
belong in local bookkeeping such as the excluded `lab-library-bindings-v1` sidecar;
they must not become scientific dependencies or be silently published. The
present merge plus immutable remote import/save CAS is groundwork for a git-like
workflow, not complete synchronization. Remote ancestor negotiation, folder sync,
deletions and divergent database heads remain separate contracts. No automatic
conflict retry or last-writer overwrite is introduced here.

## Focused fixtures

`crates/library-sharing/tests/import.rs` covers:

- A real two-run PCR/`FromRun` publication imported beside a newer unpublished
  source branch and an unrelated private plasmid, followed by fresh exact native
  material reconstruction; repeat import is idempotent.
- Strict selected conflicts and explicit checkout with a serialized complete old
  branch bookmark, plus rejection after the selected local pin changes.
- Exact-event fast-forward proof and rejection of an unrelated event with the
  same molecule revision.
- A valid same-byte attachment with conflicting metadata, a broader nonclosure
  workspace, and corrupted native revision content; failure leaves inputs intact.

No browser or HTTP checks are needed for this pure merge increment. Source is
formatted and frozen before the parent's coordinated focused native checks.
