# Native authentication, storage and concurrency security review

Integrated verification: native/WASM compilation and the connected **490-case**
focused batch pass, including the relevant admission/owner/search regressions.
[Exact case list and scope](features-security-current-integration.json).
Application HTTP, browser and live-provider acceptance remain unexecuted.

Source review and assigned follow-up fixes in the current working tree,
28 November 2026. Initial review was read-only; the bounded native decoding fixes
below are now implemented in source. No compiler, test, browser, application HTTP
check or network attack was run. The review did not read runtime databases, credential files or protected
preview directories. Findings below distinguish concrete source evidence from
unconfirmed candidates; this is not a security all-clear certificate.

## Confirmed issue: initialization can bypass hosted session quota

**Priority: medium, authenticated availability/resource bounds. The integrating
parent has corrected the classifier in source; integrated checks and focused regressions pass.**

Before correction, `crates/server/src/hosted_mcp/mod.rs:506–514` recognized an
initialization as requiring serialized session admission only when `method == "initialize"`
and `params.protocolVersion != "2026-07-28"`. The resulting flag gated both the
initialization mutex (`:521`) and `max_legacy_sessions_per_workspace` (`:529–533`).

That classification differs from the pinned transport's actual behavior. In
`rmcp 3.11.0`, `src/transport/streamable_http_server/tower.rs:545–556` explicitly
routes **every initialize request** through its sessionful handshake, regardless
of the requested protocol version. The transport defaults `legacy_session_mode`
to true (`:341`), unchanged by Collagene's workspace constructor. Its sessionful
handler calls `SessionManager::create_session` (`:2316`) and negotiates an older
supported handshake version even when the requested version names 2026-07-28.

A currently authenticated member can repeatedly submit this shape, without an
MCP session header, to one of their own workspace endpoints:

```json
{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2026-07-28","capabilities":{},"clientInfo":{"name":"quota-probe","version":"1"}}}
```

Normal authorization, fixed-lab, canonical Host/Origin and media-type headers are
still required. Each accepted handshake bypasses the host's count check and can
create another retained SDK session. The client can complete the initialized
notification using the negotiated session/version. Workspace count limits and
request semaphores do not bound the accumulated protocol sessions once individual
HTTP responses finish. This is not a cross-member document access finding; the
impact is bypassing the advertised four-session quota and accumulating session
state/tasks within an authenticated workspace.

**Bounded correction:** every sessionless `initialize` must use the same mutex
and session-count admission, independent of the requested version. Modern
stateless discovery remains a distinct method. Update the current classifier
and callers directly; do not add a second transport or compatibility fallback.

**Regression evidence:** a pure classifier case for initialize with
2026-07-28, older versions and absent version; the same cases for stateless
`discover` must not acquire the initialize path. The deliberately ignored
HTTP fixture fills the configured quota with 2026-named initialize messages
and requires rejection of the next message. Actual HTTP execution remains deferred
under the repository workflow. The pure regressions pass; the authored HTTP reproductions remain deliberately unexecuted.

## Confirmed issue: native planning collections were bounded after allocation

**Priority: medium, untrusted import availability/resource bounds. Source fix
implemented; focused regressions and the integrated compiler checks pass.**

Native Candidate decoding constructed unbounded `DesignPlan.inputs` and `steps`
vectors before `Candidate::validate` and verifier limits rejected counts above
1,000. A current Postcard `Overlap` operation with empty fragments/overlaps and
false topology occupies four bytes; about 4 MiB can therefore describe a million
large enum values despite the later semantic 1,000-step cap. The valid checksummed
outer CGP frame bounds encoded bytes, not the expanded enum allocation. The same
post-allocation pattern affected planning snapshots and variable operation lists.

The current source uses the model-owned `decode::bounded_vec` visitor on plan
inputs/steps, planning documents, all variable Source/OrientedSource/TopoSource
lists, overlap indices, the complete-digest enzyme list and the canonical
catalogue-selected enzyme list. Limits reuse the existing 1,000 inputs/steps/
sources, 128 digest-enzyme and 64 catalogue-enzyme capacities. Known
oversized length hints fail before reading an element; unknown hints use an
element seed that refuses the first item beyond the limit before constructing it.
No current JSON/Postcard encoding changes, fallback decoders or lower molecular
feature/import limits are introduced.

Three cases added to the existing history/workspace transfer suite cover:

- Valid-checksum current CGP/CGD frames with hostile claimed count prefixes,
  including 1,000,000, followed by at least that many invalid varint bytes.
  Postcard suppresses its sequence length hint when the encoded tail is shorter
  than the claimed count. Supplying the tail makes the known-count test valid:
  a custom admission error must precede decoding the first malformed element.
- Current operation field prefixes for all variable dependency-list variants,
  overlap-index and enzyme limits, plus exact existing 1,000-item admission.
- Unknown-size JSON arrays ending in a malformed extra struct after the admitted
  boundary. The resource-budget error must precede parsing that extra element.

Existing candidate/planning round trips are retained and the new count-prefix
case also confirms byte-identical current round trips. These tests were authored
and statically reviewed; this worker did not execute them. The parent owns the
frame-wide model decoder changes and final native/WASM integration.

## Reviewed fix: stale DELETE is pinned to its admitted workspace

The parent had already identified and fixed this race before this review:
request A captures a workspace, another DELETE retires it, a POST creates its
replacement under the same label, and delayed A must not close that replacement.

Current `WorkspacePool::retire` checks `Arc::ptr_eq(current, expected)` while
holding the workspace-pool lock, and `serve` passes the exact `lease.workspace`.
`Inner::remove` preserves that required pin. A stale request therefore gets a
conflict instead of removing a replacement. The source fix addresses the
reported interleaving; its regression has not been run by this reviewer.

Retired workspaces remain charged to the same quota pool until HTTP leases are
released and the scientific compute permit returns. `close_compute` closes future
semaphore admission before transport cancellation. This matters because dropping
an RPC future does not stop an already running `spawn_blocking` operation.

## Candidate refuted for the normal desktop mutation path

An initial review candidate concerned `Repository::import_tree_node`: it checks
the local parent generation before fetching remote payload bytes, leaving time
for a concurrent local rename/move before the reliable import transaction.
The reliable create/import operation itself does not carry a parent-generation
field. Examining only those two functions suggested a placement race.

The complete product path has an additional guard. Normal DesktopHost mutations
call `Repository::mutate_local_reliable_with_gate` (`desktop-host/src/library.rs`
around 138–149). `local_outcome_with_id` increments the repository epoch and
atomically writes its checkpoint for every successful mutation
(`repository/src/lib.rs:430–446`). Tree import retains its previous checkpoint
and supplies it to `state_write` in the same transaction as the import
(`repository/src/folders.rs:759–821`). A concurrent normal rename/move changes
that checkpoint; the stale import's checkpoint CAS fails and rolls back the
complete effect. Independent host and sync lanes therefore do not defeat this
existing transaction guard.

Calling bare lower-level `Library::rename_folder` through an independent trusted
connection would bypass repository checkpoint maintenance, but no supported
untrusted product path doing that concurrently was established. This candidate
must not be entered as a confirmed normal-desktop vulnerability or used to
justify a broad mutation-schema refactor without further evidence. The parent
was notified of this correction promptly.

## Confirmed issue: workspace ZIP JSON amplified before typed admission

**Priority: medium, untrusted archive availability/resource bounds. Source fix
implemented; integrated compilation and focused regressions pass.**

The previous `documents/workspace.rs` JSON boundary constructed an unrestricted
`UniqueJson(serde_json::Value)` tree before decoding typed records. The exact
feature/operation path cloned that input tree and serialized another typed tree.
Standalone feature and operation chunks decoded `Vec<Feature>`/`Vec<Operation>`,
bypassing the `Molecule` and `DesignPlan` collection guards entirely. Manifest
collection limits likewise ran after materialization. Encoded ZIP entry limits
therefore did not enforce the advertised decoded collection bounds.

The current codec parses bounded typed records first. It never constructs an
input Value. Exact feature/operation checks stream the original input against
one admitted typed serialization Value; this retains field/value equality and
nested duplicate rejection without input-tree cloning. Other typed records use
a streaming duplicate-key pass. Manifest format/version probing uses an ignored-
field header visitor rather than materializing arbitrary nested JSON.

Standalone arrays enforce 100,000 features and 1,000 steps before constructing an
extra element. Features charge one enclosing native frame (one million aggregate
features/location nodes/qualifiers, envelope text allowance), and complete
per-molecule feature admission precedes canonical comparison. Manifest vectors,
revision/blob/attachment maps, document navigation stacks, history parents,
metadata references and candidate inputs use their existing capacities. Metadata
value records read their encoding and blob as typed strings even when the blob
precedes its tag, avoiding Serde adjacent-tag generic Content buffering. Borrowed
molecule and aggregate frame admission precede save-time hashing/serialization;
reconstruction also admits complete molecules before revision hashing. Current
JSON/ZIP encoding and explicit unsupported-version handling remain in place.

Added focused codec fixtures for malformed first-excess step/feature elements,
map admission before reading an excess value, metadata content typed before
its encoding tag, nested duplicates and ignored feature fields. Full ZIP save/load fixtures retain 100,000 features, 1,000
documents and 10,000 revisions, including deterministic reopen/save bytes. These
fixtures are authored but have not been compiled or run by this worker.

## Other boundaries inspected

- Member authentication checks all configured keys, rejects duplicate credentials
  and identities, and resolves current membership again after slow body admission.
  Bearer credentials are removed before SDK request extensions receive HTTP parts.
- Workspace keys bind fixed startup lab, authenticated actor and a bounded opaque
  label. Labels cannot choose filesystem paths; session managers are per scope.
- Hosted request bodies, concurrent requests, response lifetime, workspace counts
  and retired slots have explicit bounds. Full resource guarantees still depend
  on the initialization-quota correction above.
- Managed storage uses retained directory capabilities, generated component names,
  no-follow file creation and platform-specific owned-handle cleanup. Same-user
  hostile filesystem renames are explicitly outside its OS-sandbox guarantees;
  they were not mislabeled as remote path traversal.
- Static HTTP responses use a bounded immutable asset snapshot. Request paths do
  not drive filesystem reads, and managed MCP storage is rejected when it
  intersects the static asset root.
- Native repository HTTP transport pins configured origin and lab, disables
  redirects/proxy/retries, injects credentials only after frame/origin validation,
  and bounds response reads independently of Content-Length.
- Native genome staging reserves its expected bytes before writing; append offsets,
  chunk size and total length guard its later additions. No arithmetic overflow
  was established in the reviewed staging path.

The native genome store intentionally has a separate 64 GiB per-workspace managed
reservation budget, outside the small-document store. Hosted workspace-count
limits multiply that explicit capacity; the shared-library `--max-lab-bytes`
option is not an aggregate MCP scratch-space quota. This is an operator capacity
boundary, not a demonstrated quota bypass in that store.

The integrating parent owns plan updates, fixes and the eventual focused native
and WASM checks. Unreviewed code and third-party internals are not certified by
these source observations.
