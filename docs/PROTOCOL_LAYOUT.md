# Documentary protocol plate layouts

`protocol_layout` assigns explicit construct identities to numeric wells on one
96-well (8 ×12) or384-well (16 ×24) plate. Coordinates in the native request are
zero-based; display/export labels range A1–H12 or A1–P24. Each assignment has an
explicit unique ID, one well and an optional presentation label. A well has at
most one assignment. Replicates use different IDs and wells with the same exact
source pin. Empty wells stay empty; opening a record does not place constructs.

A planned assignment binds the exact batch attachment hash, construct ID and
recorded declaration revision. An unresolved row without a declaration revision
must first be resolved in the batch workflow. Loading uses the real native batch
reader; row ordinals, names and equal DNA are never identity substitutes. Planned
rows retain their historical status, including Invalid, without becoming approved
candidates or released biological products.

A released assignment binds the exact release manifest hash, construct ID, run ID
and named output document/event/revision (`output.construct`). This can be a
reference-preserving named descendant of the recorded final product, as checked
by the existing native release reader. It is never replaced by a current head,
first matching sequence, final-step ordinal or same-name plasmid. Loading a layout
is documentary validation; it neither reruns nor adopts historical science.

`preview_revision` binds the current native workspace fingerprint and literal
request. `preview` resolves actual immutable archives and returns a bounded Layout
with source labels/statuses and exact pins. `capture` adds only its canonical JSON
attachment (`protocol-layout-v1`), preserving all document heads, history, current
candidate, old attachments and recorded assessments. Its receipt validator removes
only that attachment and compares the full original workspace fingerprint. A stale
preview cannot authorize capture under a newly computed request revision.

`load` checks attachment type/hash/size, the current schema, bound rows and actual
native dependencies. Later head edits do not alter archived assignments while the
exact original events remain retained. Missing archives or substituted pins refuse.
Ordinary workspace loading retains raw documentary sidecars; it does not grant
layout validity or physical authority. Semantic layout reading uses this module's
load/validation APIs. No recursive workspace-validator hook is installed.

Limits are384 assignments,128 layout attachments,64 Unicode characters per user
label,4MiB per layout/report/export, and32MiB aggregate referenced batch/release
archive bytes charged before decoding/caching them. Source native workspace,
release, batch and ancestry limits remain in force. Bounded collection visitors
reject excessive decoded assignments/resolved rows before reserving their claimed
capacity. Explicit invalid coordinates/IDs, duplicate wells or assignment IDs and
oversize assets refuse; nothing is silently truncated or auto-renamed.

Quantities, concentrations, mixture recipes and reaction volumes are missing for
every assignment in this first scope. The report says so explicitly. Existing
protocol arithmetic can be connected later, but no zero/default amounts, standard
recipes, well capacity, supplier conditions or experimental success are inferred.
This module introduces no scientific plan operation, sample/inventory/barcode,
reagent-stock tracking, instrument command or workflow scheduler.

CSV exports preserve explicit numeric coordinates, IDs, typed source JSON,
independent target pins, recorded engine/status and missing-quantity fields. HTML
includes an escaped plate grid and exact identity table. PDF is generated entirely
in Rust with the bundled font and has plate-specific title/producer and
Layout/Context footer labels; it asserts no scientific Status. Full hashes/events
are printed in the body. Comparison documents use the same documentary renderer
with their existing Reference/Query labels. Unsupported font/shaping is escaped
with explicit export warnings; native file acknowledgement remains the UI's job.

## Explicit publication scope

Ordinary selected-plasmid publication omits plate layouts and their batch/release
sidecars. Explicit layout publication is a separate deliberate choice: a referenced
batch or release is an indivisible immutable archive and may contain unassigned
constructs/choices. `source_bindings` enumerates the full actual dependency content
of every included archive: all design choices, existing targets, recorded inputs,
release runs and all run output events, including unassigned rows. Its attachment
list includes associated batch, manifest, run and primer-table hashes. Their typed
native ancestry closure remains the sharing layer's responsibility. No row is
silently removed or same-DNA reference substituted to make a partial archive load.
The publication reviewer must show this broader scope before disclosing it.

Nine native fixture groups use actual batch computation and atomic release, not
invented summaries. They cover96/384 boundary wells, explicit replicates and empty
plates; exact planned/released pins and equal-snapshot alias rejection; native
restart; immutable heads/history; complete indivisible dependency scope; duplicate,
label, decode/byte/count bounds; stale preview and forged capture/export receipts;
CSV/HTML escaping; independent PDF text/xref readback; and invalid historical
constructs remaining unapproved. Existing PDF readback helpers are shared directly
by protocol, ligation-outcome and plate fixtures.
