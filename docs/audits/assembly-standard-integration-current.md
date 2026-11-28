# Assembly-standard integration: current source and verification status

Updated 2026-11-28 after the collected codec correction. Selected native/WASM
compilation and 2,959 distinct effective cases pass across initial and affected
corrective executions. Not every case ran on the final binary.
[Exact execution/source stages and limits](connected-design-current-integration.json).
No browser, CDP, screenshot or application HTTP check ran. This is an
implementation/evidence inventory, not a second development plan.

## Current canonical contracts

* `planning::plan::DesignPlan.standard` is required-nullable
  `Option<cloning::assembly::standards::declaration::Declaration>`.
* `VerificationReport.standard_review` is required-nullable
  `Option<DeclaredReview>`, retained separately from status-changing `CheckFinding`
  values. Source-format failure/unknown/refusal does not mutate native physical
  reconstruction, target agreement or overall cloning status.
* The declaration pins the installed profile ID/identity and explicit ordered
  original `PlanInput` roles/IDs/revisions/intervals/orientation/slots/body bounds.
  It carries no new molecule tree, generated DNA step or alternate plan.
* `standards::declaration::assess_declared(plan, &mut DeclarationBudget)` resolves
  actual admitted plan inputs and delegates to the existing native standards
  `Request`/`assess`. Native `SourceBinding::map_interval` remains the coordinate
  owner. Stored-material offsets are not inferred into reference coordinates.
* `DeclaredReviewState` is `Assessed(existing Report)` or typed `Unavailable`.
  Its own engine, plan/declaration/review pins, requested/effective limits, scope
  and charge are retained. Documentary validation is distinct from fresh replay.

`VerificationReport::validate_standard_for_plan(&DesignPlan)` is the cheap current-plan
association guard used by source-aware native result validators, transport,
MCP and web consumers. It requires the declared review when the admitted plan has
a declaration, and does not rerun the scientific checker in a GUI render. Legitimate
global admission failures are handled explicitly without rehashing refused input.

## Exact current versions

| Contract | Current source value / decision |
|---|---|
| Worker envelope | 138, shared with F20/F28 changes |
| Candidate envelope | 35 |
| Workspace plan marker | 33 |
| Verification engine | `science-39` suffix |
| Editable assembly draft payload | 24 |
| Editable draft attachment kind | `assembly-draft-v24`; the initial freeze's stale kind was corrected in the collected diagnostic batch |
| Library recovery wrapper | schema 3 / `assembly-library-history-v3`; unchanged layout, strict current embedded draft |
| Named design wrapper | schema 1 / `assembly-design-v1`; unchanged layout, strict current embedded draft |
| Workspace container | version 2; separate plan marker handles the changed plan |
| Cloning RunRecord | **schema 2**, changed `RecordedAssessment.standard_review` retention |
| Run attachment category / separate capsule wrapper | `cloning-run-v1` / schema 1, distinct unchanged wrapping contracts |
| Overlap recipe | `overlap-recipe-3` |
| Ligation recipe | `ligation-recipe-2` |
| Tailing recipe | `tailing-recipe-2` |
| Golden Gate reaction | `golden-gate-reaction-2` |
| Combinatorial selected design | `combinatorial-source-design-5` |

Current nullable fields do not use an absent-field decoder fallback. Whole-plan
serialization binds the declaration in the current plan revision. Superseded
internal encodings are rejected through their current schema gates; no migration,
old hash path or plan compatibility adapter was introduced. Unrelated molecule,
RNA, CRISPR, TOPO and enzyme catalogue formats were not blanket-bumped.

## Registered source paths

| Owner | Authored integration |
|---|---|
| `core/cloning/assembly/standards/{declaration,declared_tests}.rs` | Declaration/review owner; source/profile/body/resource refusal; shared budget; native delegation and literal fixtures |
| `planning/plan`, verification source/engine/evidence/replay budget | Canonical optional field, admission, same native advisory pass, bounded documentary/current-plan binding, separate advisory result |
| `cloning/assembly/recipe`, `cloning/ligation/{recipe,tailing_recipe}`, `cloning/golden_gate/reaction` | Required-nullable recipe/request field carried into actual plan before hash/verification; plan-bound preview validation |
| `planning/combinatorial` | Exact optional declaration per `SelectedConstruct`, retained in resolved actual compilation; not shared across all source alternatives |
| `documents/transfer`, workspace/candidate codecs, saved material/run lineage | Current plan/report gates, source binding and recorded advisory retention; RunRecord 2 |
| `web/assembly/draft/assembly_standard_{model,view}`, preparation `MethodsState` | Required nullable editable settings, explicit profile bind, stable role/event/selector/preparation pins, blank slot/body fields, collapsed controls |
| Wizard/review/pool/library/source-choice models | Before-hash forwarding, stale binding guard, exact candidate declaration comparison, active/inactive construct retention and recovery/source inventory |
| Named saved-design model | Current draft field retention; exact standard source event inventory/attestation; no restored scientific review |
| `core/output/standard_review`, `verification/output`, `output/pdf` | Same fresh review rendered in HTML/PDF with native pins, scope, full findings, provenance and refusal detail |
| `web/protocol/verification_view` | Separate paged source-format advice and normal typed checks; native original-coordinate mapping for finding display |
| Protocol controller and saved-run review model | Exact actual-plan advisory response association before accepting a result |
| MCP assembly proposal/verification | Required nullable declaration preserved in proposed candidate; compact nullable advisory summary and full evidence artifact |

Direct current Rust constructors and SDK JSON callers were updated, including
explicit `None`/`null` when no declaration is intended. Proposal-only scientific
request types that do not carry a declaration must say `None`; a wizard candidate
guard prevents that absence from silently replacing its active declaration.

## Resource and historical integrity

Descriptor bounds precede local source cloning/hash/profile work. The existing
canonical plan/model/nested-origin admission remains in force. Original input
assets are measured before constructing the existing native request. The review
wrapper reserves bounded metadata overhead beyond the native report payload.

The native checker reports actual comparisons only when it returns a report.
`DeclarationBudget.charged()` and `observed_comparisons()` intentionally differ:
reserve the effective allowance before the call; settle to measured comparisons
only after native report and final wrapper validation. A native/late retention
exception retains an explicitly typed reservation and compact unavailable state,
not an invented measurement and not an absent review. Incomplete native searches,
IUPAC ambiguity and unsupported pins/source scopes never become conformity.

The host stores exact editable source/profile/body declarations, including inactive
constructs and previous recovery snapshots. Event inventories are checked against
real retained document/event/revision identities. Reopening restores declarations,
not a report-derived material token, ready flag or old format pass. Changed source
bindings remain editable/stale until explicit rebind/clear; numeric body offsets
are not silently reused for a different selector/choice. Reordering retains stable
identities and invokes the current grammar on the actual declared order.

## Meaningful source fixtures and check status

Authored native fixtures include RFC10 formatted additions versus body-only
profiles; RFC1000 partial role paths; both complete MoClo paths and invalid order;
reverse/circular literal BsaI findings/native source mapping; source/profile/body
and unsupported coordinate refusal; shared incomplete searches; real 4 KiB native
metadata refusal with retained reservation; and documentary/current association.

Five registered wizard workflow fixtures cover hashed-plan forwarding, native
advisory failure without product/status changes, actual stable-role reorder,
genuine same-DNA/new-event and selector/preparation staleness, explicit blank body
fields, saved-design/current-codec persistence without review authority, and fresh
HTML/PDF identities. Native/MCP source fixtures cover Some/None retention, binary
and actual SDK/artifact boundaries, unsupported versions and missing required
current fields. Their passing execution is recorded in the shared evidence, with original failures retained.

Initial live compilation identified three missed native recipe constructors in
`cloning/assembly/{edited_primers,junction_edit,selective_edit}.rs`. Those
proposal-only request types have no declaration input; their recipe literals now
explicitly declare `None`. The same correction batch fixes the attachment kind,
three web recipe fixtures, a desktop plan fixture, two selected-construct fixtures
and the screening owner's distinct fresh-source types. The initial compiler corrections pass; the later actual archive-codec omission
and its affected native/SDK cases also pass. Source formatting alone proves no
scientific result.

## Deliberate limits

The installed profiles are scoped source-part format/grammar facts, not physical
or full-standard acceptance. No mandatory cloning-status gate, automatic slot/body
classifier, source DNA/label matching, PCR-tail/clamp insertion/stripping, actual
prepared-strand/reference conversion, destination-vector or general frame/function
validator was added. Source/record/report aggregate caps can refuse large combined
editable histories explicitly; no event or evidence is silently truncated.

## Correction discovered by the connected SDK fixture

The initial caller sweep incorrectly treated `documents/workspace.rs`'s archive
loader as a new-plan producer. Its private `PlanRecord` encoder omitted the new
standard declaration, and its decoder set `DesignPlan.standard` to None. A saved
Some declaration therefore failed candidate revision validation on MCP archive
import, before session restoration. The earlier source-review coverage claim was
incorrect: updating public Plan serde and StoredCandidate did not cover this
separate current archive codec.

The current correction adds a required-nullable declaration to PlanRecord,
saves the actual validated declaration and moves it into the decoded Plan after
borrowed bounded admission. The unreleased plan schema remains33. A focused
archive fixture checks Some and None retention and hash-valid missing-field
refusal; the existing real SDK fixture covers archive/session/candidate export.
These corrections pass native/WASM checks. Five affected archive cases and nine
SDK cases pass, including the new current-format archive fixture and full
Some-declaration SDK round trip. The prior six failed core fixtures also pass;
the exact evidence retains their original failures and source stages.
