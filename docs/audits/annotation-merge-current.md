# Non-destructive annotation merge: current implementation

28 November 2026. The F19 core, ReferenceEdit, worker, web and MCP integration passes native/WASM checks, 500 distinct effective cases and the default memory-only real MCP stdio check (230 tools). [Exact sources, original failures and execution limits](annotation-merge-current-integration.json). Not every case ran on the final test binary. No GUI or provider acceptance is claimed.

## One current decision owner

[sequence::annotations::merge](../../crates/core/src/sequence/annotations/merge/mod.rs) owns explicit conflict review, contribution storage and atomic application. Existing library scanning and GFF3/GTF/BED import/export remain the source owners. Superseded `library::apply_selected*`, `formats::annotations::apply_annotations*` and ReferenceEdit `ApplyFeatureLibrary`/`ImportAnnotations` implementations are deleted. There is no hidden-Append wrapper or second duplicate policy.

| Public contract | Behavior |
| --- | --- |
| `Request` | Schema 1, exact target revision, one `Source`, 1–1,000 distinct `Selection { proposal_id, decision }` values and explicit limits. Source/input preflight precedes large copies and recursive serialization. |
| `Source::Library` | Exact library and source-bound scan report; the existing scan kernel reconstructs evidence once for the whole merge. |
| `Source::Dataset` | Exact annotation import report; the existing source bytes/hash/format/seqid owner reconstructs it once. Actual proposal IDs are used, not positional guesses. |
| `Source::Donor` | One donor molecule, declared EventBinding, explicit reverse/origin transform and selected original feature ordinals. Actual complete DNA and declared topology must match the target under that transform. |
| `revision(target,request)` | Binds the engine and complete current request, including source, decisions and target content pin. Structural validation does not scan or import again. |
| `preview` | Returns source completeness, selected actual proposed features/contributions, exact target annotation pins, typed relationships, comparison work and completeness. It creates no document or authority. |
| `apply` | Reconstructs the source once, validates every decision, computes prospective changes and publishes `Applied { molecule, review, receipt }` atomically. A late invalid decision returns no changed document. |
| Documentary validation | `Preview::validate_binding/validate_request` and `Applied::validate_binding/validate_request` check bounded current content/context/receipt and exact source/output binding without another matcher. `Preview::validate_for` is the explicit heavy fresh reassessment route. Imported previews are never application authority: apply reconstructs its own source. |

`Source::from_library`, `from_dataset` and workspace-aware `from_donor` admit borrowed data before copying it into an owned request. `proposal_ids` exposes structural source IDs without a matcher. `Source::validate_workspace` resolves a donor's actual retained event and compares the complete model; a standalone molecule/EventBinding tuple only declares a documentary association.

## Explicit decisions and relationships

`Decision::Append` adds a new annotation marked `needs_review=true`, preserving all existing feature ordinals and DNA. Exact/equivalent duplicates, including two equivalent Append decisions in one batch, refuse; the caller must explicitly choose Skip or a valid attribution decision. Mere overlap or conflicting annotations are not silently collapsed.

`Decision::Skip` records an explicit decision without changing annotations or claiming a new contribution. Even an all-Skip ReferenceEdit remains a typed history event rather than disappearing as an unrecorded no-op.

`Decision::AddProvenance { target_index, annotation_revision }` targets one exact existing annotation. It requires equivalent ordered biological geometry, kind and semantic qualifier values. It changes only the contribution qualifier: biological qualifiers, annotation order and `needs_review` remain unchanged. An already-retained identical contribution is idempotent and is recorded as not newly added. It never reconciles contradictory meanings or clears review.

The relationship owner uses canonical `BiologicalPath` order and orientation, not labels, bounding boxes or a sorted genomic interval. Semantic comparison excludes the one contribution ledger and intentional external `annotation_exchange:source_record` evidence; qualifier order/content otherwise remains exact. Review state alone is not a different biological meaning.

Typed relationships distinguish whole-feature `ExactDuplicate`, `EquivalentAnnotation`, `SameGeometryConflict`, `OppositeOrientation`, `Overlap` and `Unresolved`. Same labels at different locations remain separate. Joined paths, reversed paths and origin-crossing paths retain biological order. Fuzzy/repeated/unsupported paths remain unresolved rather than being reduced to broad intervals. Library `same_label_geometry_hint` is explicitly a coarse display hint, not precise duplicate identity or permission to apply; the old `already_annotated` field is removed.

Incomplete target relationship review refuses mutations; explicit Skip remains available. Incomplete library search is retained separately: known selected proposals do not establish absence, unique correspondence, complete source search or biological function.

## One typed contribution qualifier

The reserved qualifier is `annotation_merge:contributions`. [Ledger](../../crates/core/src/sequence/annotations/merge/provenance.rs) is its sole current parser/writer: `Ledger { schema_version:1, contributions:Vec<Contribution> }`. Unknown schemas, duplicate reserved fields, malformed content and duplicate contributions refuse. `Ledger::from_feature` provides the same bounded attribution inspection API to native, UI and MCP consumers.

A contribution retains the semantic proposed-feature content pin and a typed `SourcePin`:

- Library: exact library/entry revisions, library name/version, source molecule revision/feature ordinal, proposal ID and concrete/ambiguous/mismatch/template-length facts.
- Dataset: exact original-file hash, external format, explicit seqid, proposal ID and physical source rows.
- Donor: exact declared event, original feature ordinal/content pin and explicit transform.

Library proposal creation now emits this ledger directly; the superseded internal `library_*` qualifier collection and document-level last-application metadata are not recreated. Original external source-row evidence remains intentional file-format support. Dataset raw bytes are retained once under the existing content-addressed document metadata, not copied into every contribution. A newly projected donor annotation references its actual donor event; it does not copy the donor's previous ledger as if those nested contributors had been independently admitted again. Existing target contributions are preserved.

The ledger is documentary attribution. Parsing a valid ledger or retaining an imported source hash does not authenticate function, sample identity, event ownership or native annotation ancestry. Exact source evidence remains necessary for interpretation and any fresh scientific check.

## Stable evidence versus operation proposals

Library scan proposal IDs bind an operation to an exact target revision.
Attribution changes that revision, so an actual fresh library scan produces a
different proposal ID without supplying additional source evidence. Dataset
proposal IDs already derive from the source-file hash, format, seqid, rows and
feature content; an actual same-file reimport retains the same proposal ID.

`Contribution::same_evidence` compares `feature_revision` and all stable source
fields, excluding only Library/Dataset `proposal_id`. Library identities,
entry/source ordinal, match counts and template length remain significant;
dataset raw hash, format, seqid and source rows remain significant. Donor event,
feature and transform identity is compared in full. Exact report `PartialEq`
still compares every proposal field for fresh-replay/tamper checks. Ignoring the
redundant dataset proposal ID in evidence membership also prevents a forged
alternative ID from disguising duplicate stable source facts.

`Ledger::contains_contribution`, duplicate validation, addition and applied-report
idempotency use this same predicate. Repeated AddProvenance preserves the first
stored contribution, including its original proposal pin, and returns
`contribution_added=false`. A duplicate stable contribution is invalid even when
its proposal ID differs. Different actual entries/files/events remain separate
documentary evidence; none authenticates biological function or native ancestry.

Two independent native fixtures use actual fresh library scans and same-file
reimports after attribution, check changed library IDs and stable dataset IDs
with an unchanged single-contribution ledger, preserve the first context and biological qualifiers,
reject duplicated ledgers disguised by new proposal IDs, and exercise explicit
Skip. These checks pass in the coordinated focused execution; the dataset fixture was corrected to expect its already stable source-derived IDs.

## Direct donor coordinates and native receipts

The declared donor transform reverse-complements first when requested, then rotates by the declared origin. A linear donor permits origin zero only; circular transforms require a valid original length. Actual calls are checked directly under that declared full-record correspondence, with no matcher or uniqueness inference. Canonical model projection runs once and preserves selected original source ordinals. Selected donor features require exact supported biological paths and one complete projection; uncertain, missing or multiply emitted projection does not become guessed geometry.

The sole typed edit is `ReferenceEdit::MergeAnnotations { request }`. Its native execution emits `AnnotationIdentity` during the same mutation: old target annotations are retained, appended ones are introductions, and attribution-only changes preserve original target identity. Added donor annotations do not inherit donor DNA or annotation identity merely from equal sequence, labels or geometry. Donor event dependence is retained by `Operation::source_event` and actual workspace/history guards; a history parent is computational dependence, not physical evidence.

The native receipt records exact decisions, destination ordinals, whether a contribution was added, the complete new-feature suffix, product revision and actual cumulative comparison work. Documentary guards validate action consequences, current ledgers, candidate content pins, ordinal rosters and declared completeness. Exact request validation additionally checks unchanged DNA/topology/name, original annotation meaning/review/order, prior contributions, metadata and contribution idempotency. No output-array matching reconstructs ancestry after the fact.

## Resource and scientific limits

| Bound | Current policy |
| --- | --- |
| Request and report | At most 32 MiB each; output may be requested from 4 KiB through 32 MiB. Combined worker source/product envelopes use explicit large-message admission. |
| Decisions | 1–1,000 distinct selected proposal IDs; donor feature ordinals also bounded at 1,000. |
| Contributions | At most 32 per feature and 64 KiB of typed ledger content. Shared raw source bytes are not embedded per feature. |
| Comparison work | Generic API default 1 million; web review uses the existing 200-million shared ceiling so accepted library scans leave room for relationship review. Hard maximum remains 200 million. Actual library base comparisons or donor call checks enter the same merge ledger as ordered relationship comparisons; application receipts include subsequent intra-batch duplicate checks. |
| Retained relationships | Default 10,000, hard 100,000 across the request. Omitted relationships are explicitly incomplete, not absent. |
| Existing source owners | Library 8 MiB/1,000 entries/10,000 hits and native 200-million-comparison maximum; external annotation input 8 MiB/10,000 rows and current projection/native tree/string limits remain. |

A known library report requiring more comparisons than admitted refuses before fresh scanning. A forged low count cannot bypass the cap: the same existing scan kernel checks the enclosing actual allowance while retaining original declared options/proposal identities, then fresh equality rejects inconsistent evidence. Work is not reset per selected proposal. Target mutations are preceded by borrowed prospective-output admission including the target, changes, one review, receipt actions and shared dataset source; native save/report guards validate the final result. Bounded decoding applies to new selection, candidate, relationship, receipt, contribution and donor-index collections.

This scope excludes generic qualifier reconciliation, partial-homology/indel donor transfer, inferred unique biological correspondence, functional validation, review clearance and physical material authority. No probabilities or external scientific programs are introduced.

## Integration and verified coverage

Worker schema 140 replaces the superseded application jobs with `PreviewAnnotationMerge` and `ApplyAnnotationMerge`; immutable native MCP operations and the existing library/dataset panels use the same explicit choices. Direct donor UI/MCP resolves an actual workspace event. Event/source/settings/epoch guards prevent stale preview/commit/navigation. The attribution view calls the typed ledger owner rather than parsing ad-hoc qualifier strings.

Native machine-JSON INSDC handling validates a ledger before export and treats physical qualifier wrapping as wrapping, preserving escaped JSON string contents without inventing newline characters. Intentional GFF/GTF/BED source-record evidence excludes the documentary merge ledger from its biological qualifier hash, so attribution-only additions do not falsify external row evidence. Text-export warnings distinguish retained documentary qualifiers from complete source/workspace/search evidence.

Passing fixtures cover atomic mixed Append/Skip, late invalid target pins, existing/independent donor contributions, actual library-rescan and same-file-reimport idempotency, preserved review state, precise overlap/orientation/geometry, literal reverse/rotation/join mapping, incomplete relationship budgets, source/report tampering, native annotation identity, workspace/history/export roundtrip, and actual worker/MCP routes. Existing scientific library/file tests use explicit current decisions; no superseded API-only test is retained as a compatibility requirement. The execution record distinguishes earlier passing binaries from the affected-family rerun and final one-case fixture correction.
