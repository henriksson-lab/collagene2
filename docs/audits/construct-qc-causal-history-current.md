# F36: causal junction QC against retained history

2026-11-28 source audit. No compiler, tests, browser or application checks. This is
an implementation contract for the next coherent batch, not a claim that causal
junction QC is already implemented.

## What can be reused now

| Current owner | Real evidence available | Limit |
| --- | --- | --- |
| `material/lineage.rs::RunRecord`, `load_run`, `replay_run` | One current `DesignPlan`, explicit input roles, independent target, original step ids and exact input/product `EventBinding`s | Historical assessment is not current authority. Resolve the exact selected event and freshly reconstruct it. |
| `verification/engine/{execution,operation}.rs` | One fresh execution owns indexed products and private material; selected-path physical operations and actual PCR/digest geometry | The actual coordinate mappings are consumed locally; no source-indexed native mapping receipt is retained for causal QC. |
| `verification/projection.rs` and `model::annotations::MappingSegment` | Real source-to-product coordinate segments passed at native construction time | `Projection` returns feature values but loses original feature ordinals. Equal labels/DNA cannot restore correspondence. |
| `cloning/assembly/junctions.rs::inspect_plan` | Fresh selected final junction boundaries, windows, native method geometry, source-specific Golden Gate paths and CDS windows | Final-step inspection, not general parent-to-product causal attribution. Existing scope and output completeness flags must survive. |
| `assembly/intent.rs::PlacedFragment`, `edit_intent.rs::AppliedEdit` | Exact original oriented contributions, source ids/revisions, retained source spans and explicit indel geometry when rebuilt by their owners | Serialized target metadata is not a retained-event receipt. Fresh build can establish declared counterfactual reference intent, not prove an actual edit-history operation. |
| `model::Molecule::replace`, `annotations/editing.rs::Edit` | Canonical edit application and endpoint-envelope annotation updates, including review/deactivation | `Edit` and retained-segment mapping are private; current history commits do not retain the typed replacement. |
| `model::genbank_edit::VerifiedLog` | Genuine independent replay of typed edits against a raw GenBank baseline/current sequence, with bounded citation correspondence | Citation-scoped, format-dependent and not bound to all native before/after event identities. Do not promote its metadata key to general edit-event authority. |
| `documents/history.rs`, `workspace.rs::RevisionEvent` | Immutable event id, document id, content revision and ordered source parents | `operation` is descriptive text. Current generic `commit(next, label)` cannot prove which sequence edit was applied. |
| Web `viewer/sequence/sequence_edit_transition.rs` | Actual local before/after pins for caret continuity | Ephemeral UI navigation receipt lacks replacement DNA and is not persisted scientific history. |

`construct_qc::JunctionClaim` currently supplies a declared numeric boundary and
optional documentary history ids. Its current-sequence register finding must not
be relabelled as causal. Causal work needs the missing native correspondence and
typed history transition below.

## One current history contract

Introduce a **required typed transition** on the current native event, retaining
the operation string only as a presentation label. Bump the native workspace/event
encoding and reject superseded native encodings; no optional-default field,
legacy-event branch, label parser or migration.

Move the current MCP `tools/sequence/reference_edit.rs::Operation` into the core
reference-edit owner, with its `apply`/preflight behavior calling existing model and
sequence functions. All editor/MCP/native callers use that one type directly.
At minimum it already describes replacement, origin rotation, reverse complement,
extraction, annotation edits, rename and topology change. Batch curation must use
its actual existing annotation owner, not an invented full-document diff.

A concrete event transition contract is:

- `ImportedReference`: an explicitly imported/opened independent reference;
  supplies no causal edit assertion.
- `ReferenceEdit { source: EventBinding, operation: core_reference_edit::Operation }`:
  canonical application to the exact immutable before event must reproduce the
  after molecule, including annotation/review/provenance policy.
- `CloningStep { run_id, step_index }`: resolve the current retained run record,
  exact original input roles and selected product event; freshly replay its plan.
- `Fork { source: EventBinding }`: exact reference-preserving fork with a distinct
  document/event identity.
- Existing derived-analysis producers reference their actual current typed
  retained record and source events, using an explicit producer kind. Unsupported
  causal interpretation remains an explicit scope outcome. This is a real import
  or derived-analysis provenance category, never an `UnknownLegacy` fallback for
  existing edit callers that have not been converted.

The event hash includes the transition. No sequence edit is inferred by comparing
snapshots, choosing equal-content ancestors or parsing labels. Native edit commit
must execute the typed operation, create the event and retain its source binding
atomically. Source/head/epoch guards, no-op policy, undo/redo branches and current
workspace transaction behavior remain owned by their existing controllers.

Actual current edit entry points requiring conversion include
`WorkspaceSession::commit_document`, `viewer/sequence/editing.rs`, clipboard
replace/cut/paste, sequence modal editors, annotation controller/feature tools,
primer annotation storage and MCP `reference_edit`. Current callbacks that accept
only `(revision, edited Molecule)` must carry the actual operation at the point
where it is known; do not recover it afterward from the edited bytes.

## Fresh correspondence receipt, not a second engine

Add a private `verification::ancestry` owner consumed by the existing `Execution`.
Suggested native types are `ReconstructionTrace`, `StepTrace`, `SourceProjection`
and `FeatureProjection`. They are not deserializable scientific authority.

Each source projection retains the actual `Source` handle (input id or original
step id), source revision, occurrence/copy id, mapping segments, destination step
and explicit interpretation: retained correspondence, template-directed synthesis,
introduced literal/substitution or homology collapse. Sequence correspondence is
not atom/strand inheritance; existing native tracked owners retain that separate
chemical responsibility.

Extend the existing `model::annotations::Projection` with source-feature ordinals
aligned to the features it actually emits. `verification::projection` combines
those ordinals with the actual destination append indices and `Source` handle.
This preserves identity when two inputs have identical DNA/revisions or labels,
and when one source contributes multiple copies. Do not add feature-label or
sequence-similarity matching.

Capture the real mapping where the current operation already constructs it:
PCR/native template orientation and intentional substitutions; digest fragment
selection; end treatment/repair/resection; tailed/annealed strands; ordered
ligation; overlap/homology collapse and Golden Gate selected fragment paths.
Reference-only/recombination/TOPO domains need their genuine existing correspondence
owners, with explicit unsupported or ambiguous cases. Missing mapping cannot be
filled from nominal cut marks, displayed gene arrows or a newly fitted alignment.

The fresh trace and typed findings come from **one native execution**. Failed or
blocked steps cannot leave a successful projection/product receipt; useful
withheld evidence remains scoped as such. Do not call the verifier once per CDS
or reconstruct the same assembly a second time to manufacture tracing authority.
Refactor `assembly::junctions` to consume that same fresh receipt where necessary,
keeping `inspect_plan` as the genuine standalone convenience over the same owner.

For biological coding paths, expose the existing exact range/complement/join atom
traversal and ordinal mapping from the model annotation owner as a bounded current
API. It already supports biological ordering for coordinate-bearing qualifiers.
Reuse it in QC; do not add another recursive location interpreter or codon decoder.

## Causal QC request and verdicts

A concrete `construct_qc::causal_junctions` request selects the exact current
`EventBinding`, an explicit retained run/transition origin and bounded parent/CDS
correspondence declarations. Load retained current records, resolve immutable
source events, replay canonical operations and compare the complete resulting
molecule with the selected current event. A run selection names its original
step, not the nth successful product. For an intermediate retained step, inspect
its actual operation, not a later final junction or a truncated replacement plan.

`CdsCorrespondence` pins the parent event and feature ordinal/annotation revision,
the exact destination feature/path claim, genetic code and rationale. Native
feature correspondence must substantiate the relationship. These are explicit
role assumptions, not verified biological annotation. Keep current source ids,
occurrence ids and exact before/after coding coordinates in each result.

Replace the misleading documentary history-id fields on the one current
`JunctionClaim` with an explicit current basis: either a declared boundary/register
predicate or a retained typed causal origin. They are different checks, not two
parallel plan encodings. Update current options/report versions and callers
straight away; keep no old-field adapters.

The causal owner reports independent facts and uncertainty:

- **Frame change:** count inserted/removed calls along the mapped biological CDS
  path and track the downstream register. Intronic/noncoding edits do not become
  coding indels merely because the genomic coordinate lies inside an envelope.
- **New stop:** compare actual codon call identities and grouping under the same
  canonical genetic-code model. A stop already present in the corresponding
  parent codon stays pre-existing even if its absolute coordinate moved.
- **Frame-dependent downstream stop:** attribute a new grouping to the retained
  upstream transition(s), not the closest junction. A +1 then −1 edit can restore
  phase downstream; preserve both causes where one codon depends on both.
- **Fusion across independent parents:** require explicit parent coding-register
  and destination-path assumptions. A single ancestral CDS may not exist; do not
  invent a before protein or a unique causal deletion from an arbitrary join.

Retain all independent blocking findings. Repeated/split/partial/fuzzy/ordered
paths, multiple occurrences without unique declared correspondence, lost source
bases, ambiguous codons, exceptions and changed genetic-code context remain
indeterminate where they prevent a specific predicate. A projected annotation's
review flag is not cleared by replay. Exact sequence-change facts may be retained
as informational evidence while annotation-role-dependent verdicts remain
unresolved. No expression, fusion function, experimental repair or unique mixture
claim follows.

The audit report must distinguish: actual recorded native transition replay;
declared counterfactual assembly-edit intent replay; unsupported/unretained
history. The latter two cannot be presented as an actual recorded edit cause.

## Batch partition, schemas and fixtures

1. **Native correspondence/model owner:** `model/annotations/{mod,editing,qualifiers}`
   plus a cohesive exact-path owner; `verification/{projection,engine/operation,
   engine/execution}` and new private `ancestry`; `assembly/junctions`. Capture the
   actual native maps and feature ordinals. Freeze the receipt API first.
2. **History/reference-edit owner:** core `sequence/editing`, `documents/{history,
   workspace}` and current event producers; relocate MCP's current operation type
   directly. Required transition/event hashing and current workspace version;
   preserved immutable source graphs and current native record validation.
3. **Causal QC owner:** `construct_qc/{mod,geometry,new causal_junctions/...}`,
   explicit current claim/report representation, fresh selected-origin resolution,
   canonical coding-path comparison and independent literal fixtures. Do not emit
   a placeholder causal pass while (1) or (2) is incomplete.
4. **Host integration owner:** shared worker jobs/results/request/response/envelope,
   web sequence/editor/clipboard/session and QC panels, MCP sequence/QC tools,
   desktop/server source handles, archive/current-schema fixtures and documentation.
   All current producers/callers change in the same source batch.

`RunRecord`/`DesignPlan` remain the one current cloning route. A schema bump is
needed wherever the changed event/claim is encoded; do not bump unrelated route
versions without an actual representation change. Restored transition records
receive bounded documentary validation; causal authority still requires fresh
replay. Captured subgraphs must use the current workspace/source-event owners,
not a new ad hoc plan or unbounded copy of the whole laboratory database.

Admission should cap selected origin and replay scope before clones (one selected
run plus a bounded explicit edit chain initially), retain existing nested-origin
limits and shared native search/material budgets, and add aggregate trace mapping,
feature ordinal and coding-path limits. Suggested interactive profile: 64 explicit
transitions, 256 CDS declarations, 1M extracted coding calls, 4,096 mapping spans,
2 MiB compact QC report. Budget/omitted mapping evidence is explicit indeterminate,
never a claim of no junction or no new stop.

Independent fixtures must cover +1/+3 insertion; deletion; pre-existing versus
new stop; an upstream frame-dependent stop far from the junction; +1/−1 register
restoration; reverse and origin-spanning joined CDS; intronic edit outside the
coding path; partial/split/ambiguous mapping; same-label/same-DNA different input
roles/events; duplicated source copies; rejected/blocked sparse original step
indices; wrong saved branch/record/event; altered target or annotations; truncated
trace; malformed record/budget admission; undo/redo/archive/reopen; and keyboard,
clipboard and MCP edits that actually retain the same canonical typed operation.
