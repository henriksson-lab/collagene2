# Ordered physical processing after TOPO transfer

`topo::processing` models explicit, ordered actions over one freshly reconstructed
transfer graph. The first actions are anchored 3′-flap cleavage and templated nick
sealing. `core::topo_actions` carries original source recipes through processing,
source-seeded contour proof and the faithful tracked-duplex bridge. Processing
accepts declarations; imported reports cannot supply chemical authority.

## Declared chemistry and conservation

`IdealDuplexAnchoredThreePrimeHydrolysisV1` is an explicit ideal hydrolysis
hypothesis. It is not a named nuclease, kit protocol, cellular-repair model or
experimental-success prediction. A cleavage declaration identifies the original
source-call bond, recorded displacement step and caller assumption. The engine
requires the actual recorded duplex-anchored flap and cuts its existing bond,
creating an upstream 3′OH and downstream 5′ phosphate. Unsupported cuts leave
the incoming graph intact and report the reason. No target is used to select a cut.

Splitting a source strand introduces internal terminal chemistry without changing
its original sequence or losing its displaced portion. A per-call partition audit
checks that every original call occurs exactly once across all product components.
Subsequent nick sealing reads actual exposed ends, pairing and template continuity.
Action order matters; a later cut does not retroactively make an earlier seal work.
All components, unpaired material, unsupported actions and remaining nicks remain
in the complete report, including when no duplex can be supplied downstream.

The independent directional fixture starts with 52 source calls. Cutting its
recorded flap and sealing two explicitly selected nicks leaves a 24 bp circular
duplex (48 calls) and a separate four-call fragment. Another literal fixture
reconstructs 54 calls as a 25 bp duplex plus four-call fragment, then digests the
actual EcoRI site and checks the frozen target and AATT ends. These are model
fixtures, not experimental validation.

## Desktop and MCP workflow

The existing TOPO review opens an ordered processing panel with its exact source
frame. A native anchored-flap proposal must be explicitly added; there is no
automatically selected cleavage or seal. Inspection exposes resulting components,
source calls, ends and nicks before the user selects a component and source seed
for contour proof. Actions can be reordered or removed. Full JSON retains the
complete graph; visible long lists are paged.

A current faithful result offers an explicit declaration-only candidate handoff.
Source/target events, orientation, request revision, workspace epoch and disposal
are guarded. Verify reconstructs the original inputs and all actions, supplying
private tracked material to subsequent operations. Step-based upstream preparation
remains unsupported until tracked modification ancestry is mapped. Portable
sequence/reference summaries cannot substitute for that state.

MCP `inspect_topo_actions` and `project_topo_actions` use the same native path.
The latter creates a restrictive annotated reference document and a full evidence
artifact. `propose_plan`/`verify_plan` support the same operation; reference opening
alone establishes no physical acceptance. HTML/PDF include ordered preparation,
cut bond/profile/assumption, seal selections, component/source seed, limitations
and readable evidence summaries. Full graph arrays stay in the JSON evidence.

## Durable evidence and compatibility

Operation20 `TopoActions` requires worker66, candidate14, workspace plan12 and
science17. Jobs150/151 process/inspect actions; results142/143 return typed product
and inspection reports. The ordinary bounded binary transport remains Postcard.
Old operation19 and its processed evidence retain their layouts and version floors.
A generic `physical_processing_reports` vector with a typed evidence enum keeps
future physical operations from requiring another report field for each chemistry.

Explicit positional migrations now cover all nine report-bearing outer DTOs
through worker65, preserving nonempty old processed evidence and inserting an
empty generic vector. Earlier pre27 bare-report migrations still run first. Older
headers reject new selectors/evidence/operations; JSON defaults are not binary
migration. RunRecord and RunReceipt layouts remain unchanged.

Retained runs atomically add `physical-processing-run-evidence-v1` attachments,
bound to the owning run, original plan and exact source/target/product events.
Missing, duplicate, orphaned or transplanted evidence is rejected. Historical
inspection preserves older inner engine labels and verifies documentary bindings,
not current chemistry. Fresh replay alone reconstructs current private material.

Worker149 and MCP historical inspection return exact old v1 bytes for old-only
runs. Runs with generic evidence return `physical-processing-run-bundle-v1`,
containing each attachment kind, hash and original JSON evidence. This also covers
mixed old/new operations without reserializing their stored graphs. Downloads
remain bounded; oversized collections use the existing whole-workspace archive.

## Validation and remaining scope

90 focused native tests pass across new processing/projection/core/verification,
old physical/nick/contour behavior, nine-path migration, retained history and real
MCP SDK calls. The application/worker WASM compiler check passes. Browser checks
remain deferred. Public port8082 still serves the previous worker65 snapshot.

Broader upstream tracked-material preparation, durable graph-only candidate review,
additional chemical profiles, other cloning methods, batch workflows and remaining
platform/delivery/GUI acceptance stay in the active project plan.
