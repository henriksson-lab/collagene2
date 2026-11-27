# Selected-plasmid publication closure

`collagene-library-sharing` derives a native, per-document publication view from
an exact `{document_id, revision, head_event}` selection. The source workspace is
never changed. It does not flatten the selected plasmid into `.cg2`, interpret a
sharing enum as permission, or fall back to uploading the entire workspace.

The crate depends only on scientific core and portable Rust serialization. It
has no database, authentication, filesystem, HTTP or external executable
dependency and can be used by native/Tauri and WASM hosts. It computes documentary
scope; existing scientific verification and physical replay remain unchanged.

## API

- `inspect_document_closure(&Workspace, &Selection)` returns an independently
  derived `Inventory`, without cloning unrelated DNA/attachment payloads.
- `project_document_closure` returns `Projection { workspace, inventory }` for
  browser/native hosts that encode through the existing workspace worker.
- `document_closure` is a synchronous convenience returning
  `Export { workspace, archive, inventory }` with the existing native archive.

An inventory contains deterministic sets of document, revision, event, candidate
and attachment identities, exact projected document heads, and omitted counts.
Omitted counts describe the source view and are not authorization evidence. The
server independently derives the permitted sets/heads from received native data
and compares them with the actual archive before granting `PublicationScope`.
It must not echo a caller-supplied inventory as proof of a safe closure.

## Supported dependency contract

The source first passes ordinary `Workspace::validate`, including molecule hashes,
event DAGs, native bounds, run schemas, exact event/source pins, nested origin
ownership and associated evidence validation. This is structural validation, not
scientific replay.

Starting at the exact selected head/event, the closure follows unchanged event
parent IDs and retains each event's exact molecule revision and document identity.
If an event belongs to a retained cloning run, that **complete immutable run** is
an atomic dependency. All its target/input/product events, original source
snapshots and reference/physical/processed evidence attachments are retained.
Pruning unused run steps would change its identity and violate native validation;
the helper does not rewrite it into a different reaction.

`ReconstructRunProduct` operations additionally name original runs and exact
selected/target/input/product event pins. Those are followed transitively. The
existing native validator proves that each origin capsule agrees with its owning
run, including the serialized candidate and nested origins. Capsule bytes are
preserved unchanged; the helper does not independently guess their chemistry or
manufacture tracked state from product letters.

Supported attachment families are `cloning-run-v1`,
`reference-design-run-evidence-v1`, `physical-processing-run-evidence-v1` and
`processed-topo-run-evidence-v1`. Companions are kept only for retained owning runs.
An event using the reserved native `Cloning run ...` producer label without its
owning record is rejected as unresolved.

Standalone workspace candidate proposals have no required foreign key from an
immutable run: the run already embeds its complete plan, and origin capsules
already contain their candidate packages. Unrelated/unpublished candidate
proposals are excluded. Ordinary molecule metadata remains exact inline content;
it grants no dependency or material authority. The helper does not chase arbitrary
hash-looking strings or treat free-form metadata as executable reference schemas.

## Local branches and presentation state

Only required event ancestry is published. Sibling/redo branches, unrelated
documents, later private edits to a source and the local workspace title are not
implicitly included. Dependency document handles are projected onto a
deterministically chosen retained event (lexical event-ID order), with navigation
stacks cleared. The selected document keeps its exact selected head/event. All
retained immutable events and molecule revisions remain byte-for-byte unchanged;
projecting a handle does not create a new event or grant new ancestry.

The publication view is named `Plasmid publication`. Its active document is the
selected plasmid. Existing native validation and archive encoding must accept
the result. The original workspace retains every branch, navigation stack and
unpublished draft.

Known workspace presentation sidecars are deliberately omitted: output visibility,
oligo-table exports, enzyme sets, named digest setups, assembly draft/recovery
records and local remote bookkeeping. The exact App contracts audited are:

| Attachment kinds | Role and source evidence |
| --- | --- |
| `assembly-draft-v14` | Current editable assembly declaration; no scientific status or physical authority is restored. |
| `assembly-design-v1` | `assembly_saved_designs_model.rs` stores named **editable declarations**; no candidate/result authority is restored. Their source pins/capsules belong to unpublished future designs, not dependencies of an existing plasmid event. |
| `assembly-library-history-v1` | Local draft-edit recovery bytes; not an immutable product producer. |
| `primer-pair-draft-v2` | `primer_pair_storage.rs` explicitly describes unfinished editable text, not successful PCR or saved binding. |
| `annotation-display-v1` | `annotation_visibility.rs` stores hidden feature indices as display preferences, without sequence edits. |
| `lab-library-bindings-v1`, `lab-library-branches-v1` | `library_bindings.rs`/`library_documents.rs` store per-remote links and retained local branch bookmarks, outside scientific ancestry. |

No generic draft/version prefix is accepted: for example, `assembly-draft-v15`
requires a separate role/schema review. Scientific feature/alignment libraries,
PCR preview/reference evidence and read evidence are not classified as private UI
state merely because the App displays them. They retain the unsupported-contract
boundary unless a dedicated scientific closure handler exists.

Recipes and oligos needed for a retained reaction remain in its immutable run plan.
Local remote bindings carry no scientific graph references or credentials and
must never be published as part of a plasmid's history.

Other attachment kinds are an explicit unsupported reference-contract boundary,
even if they might be unrelated. There is no generic attachment foreign-key
schema proving that dropping them is safe. The current helper rejects them
instead of guessing ownership, silently losing ancestry or publishing every
attachment. Read-assembly/validation, arbitrary extension assets and other
attachment contracts can receive dedicated closure handlers later. An explicit
whole-workspace publication remains a separate authenticated host choice, never
an automatic recovery from this error.

Core workspace limits apply; the index additionally caps reference traversals at
one million. Sets and queues ensure each retained event/run is expanded once.
No PCR, digest, assembly or physical material replay runs while deriving scope.
Reopening and downstream use still execute those existing algorithms against the
unchanged retained source/event declarations.

## Offline local library and later sync

This export is a publication view, not a replacement for the local library. A
Tauri/offline library must keep its complete local plasmids, folders and unpublished
branches. Git-like synchronization with selectable remote databases is a separate
typed operation; it must retain remote bindings independently for each remote,
distinguish library revision IDs from native document/event/revision pins, and
report divergent heads explicitly. No external `git` process is required or
introduced here. The current shared-library wire is not silently reinterpreted
as a complete merge or bidirectional sync protocol.

## Validation

Focused fixtures construct actual native PCR runs, reuse saved material through
an explicit origin capsule in a second recorded PCR, publish the selected product,
load its native archive and freshly reconstruct/reuse it again. Independent literal
DNA is checked after replay. Unrelated private documents, later source edits,
candidate proposals and local drafts must be absent, while exact ancestor source
events remain available. The saved-material fixture uses the actual current App
`assembly-draft-v14` attachment kind with private working bytes. A separate fixture
checks all ten known draft kinds, named designs, primer drafts, visibility and
remote/branch bookkeeping, including opaque unfinished/recovery bytes, and proves
that unrelated private DNA does not enter the reopened publication. Future draft
versions and scientific reference schemas are explicitly rejected. Separate
fixtures cover a plain plasmid, stale selection pins and unknown schemas.

The registered crate passes its focused publication and import fixtures in the
unified native graph. The native server and browser integration consume these
helpers; complete offline/multi-remote synchronization remains separate work.
See the current verification checkpoint in [the plan](PROJECT_PLAN.md).
