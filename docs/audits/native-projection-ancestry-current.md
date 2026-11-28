# Current native projection ancestry receipts

Status: integrated model/verifier/edit/causal owners pass selected native/WASM checks and focused native fixtures. [Combined compiler and execution evidence](causal-history-current-integration.json).
This is the reference-coordinate receipt foundation for causal CDS checks. It
is not an executable plan, an imported authority object, a chemical strand
history, or a biological certificate. The public causal-QC rule layer is documented in [its current audit](construct-qc-causal-current.md).

## Frozen private core API

Owner: `verification/ancestry.rs`, with independent fixtures in
`verification/ancestry/tests.rs`. All receipt types are `pub(crate)` and none
implements `Deserialize` or `Serialize`.

- `TraceBuilder::new(planned_steps, Limits) -> Result<Self, Error>`;
  `register_input(id, &Molecule) -> Result<SnapshotToken, Error>`;
  `begin_step(original_index) -> Result<StepBuilder, Error>`; `finish()`.
- `StepBuilder::destination(&Molecule) -> Result<SnapshotToken, Error>` creates
  an actual native draft snapshot. `project(ProjectionCapture)` records the
  canonical annotation projection and exact reference spans actually used.
  `introduce(token, start, end, IntroductionKind)` records explicitly introduced
  calls. `seal(token, &Molecule)` binds the actual final annotations/metadata.
  `commit(final_token)` records the successful original step slot.
- `ProjectionCapture` fields: `source`, `source_molecule`, `destination`,
  `destination_molecule`, `occurrence`, `kind`, `segments`, `projection`,
  `destination_feature_indices`, `changes`. Occurrence is explicit, not derived
  from equal molecules, labels, positions, or pointer addresses.
- `ProjectionKind`: `ReferenceCorrespondence`, `TemplateSynthesis`,
  `EditedTemplate`, `HomologyCollapse`. `MappedChange` holds exact
  `source_position`, `destination_position`, source-forward `before`, and
  destination-forward `after` calls. An ordinary copy cannot conceal a change.
- `IntroductionKind`: `PrimerTail`, `DeclaredLiteral`,
  `NewlySynthesizedUnmappedCalls`. Introductions are not fabricated template
  origins; reference correspondence does not certify physical inheritance.
- `ReconstructionTrace::step(index)`, `projections(index)`,
  `snapshot_source(token)`, `snapshot_projections(token)`, `lookup_session()`.
  `StepTrace` binds `step_index`, `product_revision`, `product_length`,
  `local_coverage_complete`, plus an opaque `final_snapshot()`.
- One `LookupSession` per QC job: `at(step, position, stop_at)` and
  `feature_at(step, ordinal, stop_at)`. Both return all captured alternatives,
  explicit `complete` and `reason`; neither returns absence when work or mapping
  is incomplete. `stop_at: Option<LookupStop>` may name a current-plan source or an exact
  opaque snapshot (`PlanSource` / `Snapshot`). Snapshot selection requires a
  sealed identified source and rejects current/future unscoped plan products. Unrelated origins remain in the result rather than being filtered out.
  `remaining()` exposes the shared remaining receipt-walking work.

`OriginCall` binds original `SourceIdentity { source, revision, scope }`, source-forward
position, composed reverse orientation, actual occurrence path, and any captured
substitutions along that path. `FeatureOrigin` binds exact original feature
ordinal and annotation revision, not a feature label. Paths preserve the actual
producing step and intermediate destination token.

A token contains an explicit private trace-owner ID and an allocation serial.
The owner ID is allocated with checked `AtomicUsize`, supported by the existing
native/WASM atomic capability; exhaustion refuses a new builder. Serial numbers
never rewind. Foreign trace tokens and tokens from discarded stages cannot alias
new snapshots with the same vector index or identical sequence. Actual intermediate
sources must be sealed before use. A step can create its final draft before an
intermediate source; sealing/use order, not creation-index order, determines the
DAG. A failed staged operation drops its snapshots and retains all resource charges.

## Current canonical model receipts

The canonical `annotations::Projection` now carries these native fields,
populated while traversing each original source feature:

- `source_feature_indices: Vec<usize>` aligned with emitted `features`;
- `unmapped_source_feature_indices: Vec<usize>` aligned with `unmapped_features`;
- `dropped_source_feature_indices: Vec<usize>` for every actual dropped original,
  aligned in count with `dropped_features`.

There is no compatibility representation. Every original ordinal must have a
consistent disposition. Repeated emissions can carry the same original ordinal;
no ordinal may simultaneously be emitted and dropped/unmapped. Receipt capture
checks these invariants, actual source/product bases, source revisions, and the
actual appended projected annotation at its exact destination ordinal before seal.

## Bounds

Default/max construction caps: 1,000 planned slots, 4,096 snapshots, 4,096
projection calls, 4,096 total reference/introduction spans, 32,768 total feature
dispositions plus changed calls, 64 MiB aggregate scanned/serialized source and
annotation evidence. Individual annotation serialization is limited to 256 KiB.
Hashing is streaming; repeated annotation hashing consumes the aggregate budget
before work. Lookup default work is 20 million, configurable up to 50 million,
shared across base and feature lookups; each individual lookup also limits depth,
pending breadth, and emitted alternatives to 64. Actual node/edge/change/path-copy
work is charged. Limits produce errors or explicit incomplete lookup evidence;
partial results cannot certify absence or a unique causal origin.

## Bounded coding-path comparison integration

The current model `annotations::biological_path` owner reuses canonical atoms
with an identity projection. `BiologicalPath::from_location` accepts `PathLimits`
and exact source length; `call_at`, `ordinal_of`, and `calls` preserve biological
order with a separate inverse interval index. Fuzzy/order/unresolved paths and
repeated identical source calls refuse explicitly. Source-call and ordinal lookup
are logarithmic in admitted span count, not a full CDS scan per base.

For an explicitly selected parent `Source` and original CDS ordinal:

1. Freshly replay the actual saved run/current design through the existing
   verifier; retain this private trace only from actual reconstruction. Bind the
   final report to exact saved source event/history identities and plan revision.
2. Resolve product-feature ancestry through `feature_at`, not label matching.
   Retain partial/split/discarded annotation evidence and review status.
3. Use canonical biological paths and one shared receipt lookup session to
   relate product coding calls to exact parent coding ordinals, retaining all
   real occurrence/role alternatives and orientation. Extra introduced calls
   and omitted source spans are explicit, not fitted by sequence alignment.
4. Use the existing canonical CDS translator, actual `codon_start`, exact genetic
   code, partial flags, and biological orientation. Compare actual codon groups
   under the captured correspondence. A moved pre-existing stop is not a newly
   created stop. A downstream stop caused by changed codon grouping can be linked
   even if it is distant from the retained junction.
5. Retain multiple causal edits, compensating +1/-1 shifts, and restored register
   segments. Do not assert a unique edit cause when genuine homology alternatives,
   incomplete coding paths, or missing receipts prevent that conclusion.

Actual typed native edit transitions are retained by `documents::reference_edit` and workspace2; an operation description string or nearby retained boundary supplies no edit authority. Compound/partial CDS and homology-collapse paths must remain
explicitly unresolved where canonical mapping is insufficient. Nothing in these
receipts independently verifies an imported annotation's biological truth.

## Independent source fixtures

Passing focused source fixtures include: identical DNA/labels with distinct source roles; exact
feature ordinals; reverse-complement palindrome through a real intermediate;
discarded-token index reuse; foreign-builder tokens; missing map; explicit
introduced literal; shared lookup exhaustion and nonresetting discarded-stage
capacity; unsealed/changed sources; actual C-to-T mapped substitution; multiple
real homology-collapse origins; appended annotation tampering. Model ordinal integration and native capture callsites are included in the combined compiler/execution evidence above.


## Same-executor integration

`verification::engine::verify_with_ancestry(plan, target)` returns the ordinary
`VerificationReport` and a separate `Result<ReconstructionTrace, Error>`. It uses
the same execution/operation owners and preserves ordinary status/products.
Trace refusal never supplies causal authority and never changes an ordinary
scientific reconstruction into a rejected result. Actual unsuccessful operations
and enclosing report/evidence withholding discard the latest trace stage while
retaining resource charges. A private trace is not a serialized report input.

The shared native projection helper returns `NativeProjection` with the actual
new-feature destination indices and the actual previous-feature dedup remap.
Every previously captured destination ordinal is remapped from that retain-loop
receipt, so earlier origins survive later appends/deduplication. Canonical helper
receipts are consumed directly, without re-running source discovery or choosing
a source by DNA/annotation equality. PCR, terminal, digest, ordered overlap,
ligation and Golden Gate paths bind actual source handles and actual same-step
intermediate tokens. Graph-owning native methods now emit their actual segment receipts. Missing or
resource-refused capture refuses the private trace, while the ordinary method
result retains its established scope.


## Retained origin namespaces and shared replay budgets

`SourceHandle` distinguishes `Plan(Source)`, canonical material `SourceRole`,
`MaterialProduct`, and original `TopoRole`. Native derivations produce a one-step
receipt graph without constructing a second plan. `OriginBinding` variants retain
exact Run selection/plan/target/input event bindings, Material product/record hash/
context/source-event bindings, or TopoGraph selection/context/role-event bindings.
An `OriginOccurrence` path preserves enclosing step, actual occurrence, and shared
immutable binding. Missing graph event provenance remains `None`; it is not filled
from equal source DNA or assumed to establish event-level causality.

`TraceBudget` is a private capability with one monotonic shared ledger across all
nested current replays. Plan slots (1,000), native snapshots/projections/spans,
annotation dispositions, namespace links (32,768), and evidence scans are shared.
`retain_model` admits source copies through the canonical model guard and streaming
byte count before retention. Moving already-owned native product/root models
avoids additional copies. Namespace depth is at most eight.

The same actual saved-run/material replay returns its private trace. Import checks
the shared capability, exact native context, and selected product; it only moves
receipts and relocates opaque token references. It neither runs a plan nor
reconstructs scientific mappings. Actual saved product presentation copies retain
the exact original native molecule pin even when its saved name differs.

TOPO expansion records aligned root scopes at each actual input push. Current
caller inputs retain their actual operand index. Saved graph roots retain their
original local role, exact graph context, real source-event bindings, and actual
expansion-node occurrence. Cached saved-material replays are captured in that same
native cache execution and imported once. Each actual graph-role presentation has
its own scoped copy edge to the cached source; lookup follows that edge to deeper
original sources, preserving graph-specific occurrence paths without cloning or
replaying the cached source. Exact graph-role snapshots can also be selected as a
baseline using `LookupStop::Snapshot`.

`ReconstructionTrace::sources()` exposes bounded source identities and their
opaque tokens for exact retained-event/role selection. It does not search names or
choose a source by sequence similarity. These receipts remain a foundation; the
causal coding/frame/stop comparison and typed edit-history adapters still require
integration using the canonical biological paths and translation owners.

## Canonical reference action namespace

The reference-edit owner uses `for_reference_edits` with the exact before-event
binding and actual canonical action revision. `register_reference_source` checks
that event and the actual source revision. Native action outputs are
`ReferenceProduct { action_index }`, with the source represented as a real
`ReferenceEvent`; no plan input or material role is invented. The same actual
output token retained after successful `NativeStep::finish` supplies the next
action's source. Canonical model/edit owners must supply actual coordinate and
feature-disposition receipts; the receipt owner does not derive edits from final
DNA differences, names, feature labels, or operation prose.
