# Web and worker ownership audit

28 November 2026. Static review of the current working tree, including the pending
coordinate-qualifier batch. Production files were not edited. No compiler, tests,
browser, screenshot, application request, helper build or release verification was
run. Line references describe this snapshot and will move during implementation.

## Current limits

`python3 scripts/check_rust_layout.py` reports **1,238 Rust files in 212
directories, zero violations and zero 3,000-line warnings**. `crates/web` contains
293 Rust files including tests; `crates/worker` has one 36-line entry point.
The largest web directories contain 16 direct Rust files (`workspace` and
`cloning/methods`); `library`, `viewer/sequence`, and `cloning/preparation` contain
15. These are not exceptions or grandfathering allowances. Regroup before growth
approaches 20.

| Current file under `crates/web/src` | Lines | Decision |
| --- | ---: | --- |
| `app/mod.rs` | 2,869 | Near warning threshold; move remaining feature behavior into existing owners before more callbacks are added. |
| `cloning/preparation/cloning_methods.rs` | 2,511 | Highest priority ownership extraction; proposals, guarded evidence and pool lifecycle have identifiable boundaries. |
| `viewer/sequence/sequence_view.rs` | 2,244 | Below warning; keep document display/selection behavior cohesive rather than splitting by line count. |
| `assembly/review/assembly_review.rs` | 2,126 | Cohesive atomic staging authority, with tests; retain owner and improve natural enum methods. |
| `assembly/draft/assembly_workspace.rs` | 1,978 | Frozen UI intent and review lifecycle are a justified view/controller boundary from serializable Draft. |
| `viewer/map/renderer.rs` | 1,755 | Renderer resources and cleanup are colocated; convert natural State operations in place. |

The old flat 228-file web root, 7,273-line App, scattered Draft identity impl,
app-owned export service and giant worker dispatch no longer describe current
source. Do not reopen them as new findings. The prior audit documents are historical
evidence, and the follow-up already records completed work.

## Ranked implementation batches

### 1. Own cloning proposals, evidence and prepared-pool lifecycles

`cloning_methods.rs` declares `PlannerForm` at 38, `Suggestions` at 72,
`GoldenGateSuggestions` at 81 and `Evidence` at 89. Its `CloningMethodsPanel`
starts at 155 and continues to 1265. It handles asynchronous proposal generation,
acceptance, explicit application, scientific evidence calculation, diagnostics,
exports and finalization. Guard logic is repeated around apply at 408/604,
diagnostic navigation at 931, suggestion navigation at 1150/1193 and proposal
export at 1221. The independent input/source/target/options/workspace-epoch pins
are important; the problem is that their lifecycle is embedded in one component.

The same file owns four large request builders at 1432 (`overlap_request`), 1495
(`ligation_request`), 1638 (`primer_design_request`) and 1757
(`golden_gate_request`). They interpret `MethodsState` and the proposal form,
validate settings and construct concrete core requests. Put those operations on
the actual preparation or proposal owner. A stateless assembly of several
independent operands may remain a named free function; adding a receiver merely
to rename a helper is not the aim.

A complete pool owner is already visible at 1900–2464: `PoolContext`,
`PoolReceipt` and its methods, `PoolSettings` and its methods, `PoolEvidence`,
`PreparedPool`, cleavage details and candidate details. Move these together into
a semantic `pool/` owner, including capture/admission guards now at 1097/1109.
Do not separate the DTOs from the meaningful receipt/settings methods. Preserve
explicit selection, intended versus alternative products, cancellation,
source/options/target/epoch binding, and independent finalization checks.

Suggested allocation: use a `cloning/preparation/proposals/` family for overlap and
Golden Gate proposals with their declaration, request construction, fresh-context
checks and explicit apply behavior, and a `cloning/preparation/pool/` family for
the pool lifecycle. Keep the panel as composition. Reuse the existing form owner
when appropriate; do not create competing copies of MethodsState. The existing
`MethodEvidence` component at 1270 has a deliberate view-type erasure boundary
(comment at 1267); presentation can remain separate from guarded state ownership.
Do not replace the panel with a universal bag of signals or distribute one new
controller's impl across many small files.

### 2. Finish canonical import paths and natural draft/intent methods

Two small ownership facades remain and should be removed in the same relevant
batch, updating consumers directly:

- `cloning_methods.rs:1429` publicly reexports
  `cloning_preparations::prepared_rows`; `assembly_workspace.rs:1616` calls through
  that component module. Import from the real preparation owner, then make the
  panel's import private. Most callers already use the canonical path.
- `assembly/draft/model/mod.rs:54` reexports core `intent::MaterialDomain`.
  `restriction/digest.rs:100/103` and `digest_fragment_model.rs:27/239/781` reach
  the core declaration through the UI draft model. Import the core owner directly.

`Draft` and its main implementation are now together in
`assembly/draft/model/mod.rs` (declaration 124, impl 793). Keep them there.
The same owner still has free functions `set_digest_fragment` (1153),
`inputs_request` (1290), `edit_request` (1420), `orient_fragment` (1459), and
`set_material_origin` (1496). Review these as natural Draft operations, convert
where they express its invariants/transaction, and update callers. Their current
colocation means no physical split is required.

`IntendedReference` at `assembly_review.rs:10` owns `target` and `fragments`
methods, but `assembly_workspace.rs:61/69` still matches its variants to derive
context and warnings. Add those accessors beside the enum, then update direct
callers. The staging functions genuinely combine a workspace, candidate and source
bindings and need not all be forced onto this enum.

### 3. Reduce remaining feature behavior in App before further growth

The recent restoration refactor is present: `apply_restored` at 887 advances the
workspace epoch and delegates reset/install behavior to narrow owners in one
batch. Do not redo it or claim it remains a monolithic restore implementation.

App still defines feature operations after composing those controllers. Concrete
examples are saved-enzyme-set persistence callbacks at 1797–1846, retained digest
setup operations at 1847–1974, antibody publication/selection at 1975–2097, RNA
save/load/delete at 2098–2158, and comparison save/reopen/admission at 2159–2246.
Move each coherent callback/state dependency group into its existing restriction,
biology or comparison controller, with narrow dependencies, leaving App to compose
the result. Keep scientific event pins and workspace epoch explicit across awaits.
This is a growth-prevention priority because the file is only 131 lines below the
warning threshold; no size violation currently exists. Modal/shell composition
from 2444 onward is a legitimate App responsibility.

### 4. Improve renderer member ownership in place

`MapRenderer`, private `State`, GPU resource cleanup and input/rendering lifecycle
share `viewer/map/renderer.rs`; that is good cohesion. `advance_gesture` at 1254
mutates exactly State's drag, range preview and rotation, while `draw` at 1497
uses its resources/caches. Convert these to private State methods and review
`feature_spans` (1296) and `nearest_span` (1319) as State queries. Update the few
calls in MapRenderer directly. Retain pure geometry/tessellation/color algorithms
as free functions. There is no reason to create separate files for every primitive
or move the State declaration away from its resource lifetime implementation.

## Boundaries that should remain

- Serializable Draft and pure validation are independently usable in native model
  tests; frozen intent/review and Leptos lifecycle are separate responsibilities.
  Registered test modules under `assembly/draft/model/tests` are appropriate,
  including explicit `#[path]` discovery. They are not unregistered coverage loss.
- The 36-line worker entry is a browser message adapter over authoritative core
  `AnalysisSession`. The new core execution directory has typed stateless domain
  operations and a bounded exhaustive router. It should not become many scattered
  impls on one universal context merely to satisfy a methods preference.
- `CachedWorker` and its meaningful implementation/drop lifecycle are together in
  `platform/worker.rs:152/172/307`; no additional split is justified by this audit.
- Same short type names such as Capture, Proposal, Context, and SourcePin belong to
  distinct modules. A name-only scan cannot establish scattered implementations;
  do not merge unrelated owners based on those names.

## Integration and acceptance

Apply each coherent ownership batch across source, direct imports/callers, schemas
if actually changed, and fixtures before compiling. Preserve current serialized
bytes and scientific behavior for organization-only changes. Remove old paths;
no forwarding facades, aliases, fallback encodings or extra plan representations.

Run the repository layout check after the completed edits. If tests move, preserve
explicit Cargo registration and compare before/after executable inventories. Build
selected integration harnesses and the web unit harness together with the same
packages/features, then execute only affected cases. Run the required WASM check
in the background while independent implementation work continues. Reuse the
shared cache. No automatic release verification or browser/application checks.

The helper's read-only ownership/coupling reports may supplement this review,
using the WASM target for gated UI code. Its automated write/release/rollback
workflow conflicts with repository policy and was not used here. This audit is
source evidence, not a semantic-tool proof or graphical acceptance result.
