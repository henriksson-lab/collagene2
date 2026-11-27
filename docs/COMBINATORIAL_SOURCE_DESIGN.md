# Alternative source bins and selected construct matrix

The assembly dock now supports ordered source bins containing interchangeable
choices. Each selected matrix row is a separate intended construct, with a stable
identity and its own editable preparation and junction settings. This is distinct
from enumerating possible products of one Golden Gate reaction pool.

## Existing dock, explicit source choices

Enable alternatives from the ordinary assembly draft. Its existing slots become
ordered bins and its current design becomes the first selected construct. Choose
a bin, use the existing source/region/orientation and preparation controls, then
explicitly capture the active slot as a new immutable alternative. Original exact
source events and actual source models remain retained.

The matrix counts and pages the Cartesian choices without allocating the full
product. Pages contain20 combinations and support a decimal page jump. Add chosen
combinations explicitly to a selected construct list, displayed10 rows per page.
Repeated combinations are permitted as distinct construct IDs, so identical source
choices can have different intentional primers or junction edits. Choosing a row
opens it in the same dock; there is no second preparation wizard.

Saving or switching constructs preserves its raw primer/tail fields, method,
topology, reaction settings and explicit junction edits. It does not replace
manual primers with suggestions. Source alternatives retain their original raw
preparation defaults; new constructs begin from those choices, while saved
constructs keep their own overrides.

Working source edits may temporarily differ from the selected immutable choices.
They remain editable and persist, but saving/switching constructs or scientific
preparation is blocked until each changed source is explicitly captured. A pending
source cannot inherit the old choice's identity. Removal refuses choices referenced
by retained constructs. The bin structure is fixed when enabling alternatives in
this initial connection; ordinary slot addition/reordering/removal is disabled
while the library is attached. Reviewed bin insertion/removal/reordering and raw recovery are now implemented
([details](ASSEMBLY_LIBRARY_EDITS.md)).

## Identity and scientific reconstruction

Stable design, bin, choice and construct IDs are separate from labels and row
positions. Source PlanInput IDs bind the complete stable tuple; fragment roles
bind their bin. A separate declaration revision includes the actual choice payload,
method, exact events and independent target. Equal DNA, repeated labels or duplicate
choice vectors do not collapse construct identity or reuse a fresh receipt.

All existing preparation projections, primer/junction proposals, reaction capture
and candidate checks use the same resolved source/role table. Legacy ordinary
drafts retain their explicit ordinal defaults. Switching constructs invalidates the
intended-reference review and method evidence even when their DNA is identical;
ordinary active-slot navigation does not invalidate an unchanged intended reference.

`core::combinatorial_design` provides bounded typed designs, checked mixed-radix
count/page helpers, explicit selected rows, actual workspace-event resolution and
existing compiler/verifier dispatch. Targets are either an explicitly selected
retained reference event or a declaration from original source contributions and
literal junction edits. Targets never come from reconstructed prepared products.
Historical event resolution leaves current document heads unchanged. The web dock
uses its existing independently declared intent and worker-backed preparation path.

The core selected-source profile is checked before copying relevant models. Native
validation remains available for unrelated rich workspace documents; those are not
forced through the smaller scientific-operation profile. Incomplete unselected
preparation choices do not block a separately valid selected construct.

## Persistence and bounded state

Editable data now lives in UI-independent Rust modules for draft, preparation,
source identity, roles and alternatives. Scientific receipts remain private and
transient. The native attachment writer uses `assembly-draft-v14`, JSON version14;
the reader explicitly migrates versions1/2 and accepts version3. Attachment-kind and
payload version must agree. Multiple active known drafts are an error. Old readers
can preserve the unknown new attachment kind instead of interpreting alternatives
as an ordinary ordinal draft.

Restoration checks every retained captured event, including inactive choices,
against its actual document/event/revision and source snapshot. Autosave removes
all known active kinds and writes only version3. Per-construct edits contain
settings and choice IDs, avoiding recursive draft snapshots and source copies
for every matrix row. Pending working source edits are retained separately from
the last choice-consistent saved row.

Limits are12 bins,64 choices per bin,256 total choices,256 selected constructs,
the existing32 MiB draft budget and existing source/model limits. Checked `u128`
rank/count arithmetic supports matrices larger than a32-bit address space without
materializing them. Stable monotonic IDs are not reused on removal. Malformed
working or saved role tables fail rather than being relabeled by index.

## Remaining delivery work

This checkpoint connects explicit choice matrices to ordinary per-construct editing
and verification. It does not complete F13/F27 bulk delivery. Batch evaluation with
progress/cancellation, immutable batch records and dedicated MCP design operations
are now implemented ([details](COMBINATORIAL_BATCHES.md)). Matrix-owned per-construct
run/output records, bulk sequence and oligo/protocol release and staged product
inspection are implemented ([outputs](COMBINATORIAL_RELEASE.md)). Finite reaction
composition now has explicit per-construct discovery and key-bound selection
([reaction matrix](COMBINATORIAL_POOLS.md)). Library topology editing now has reviewed transactions and recovery
([details](ASSEMBLY_LIBRARY_EDITS.md)). Final acceptance remains in scope; extend the typed design layer and existing compilers. Fusion-end ranking and domestication
remain distinct planned capabilities.

Browser/visual/performance acceptance remains deferred until all in-scope code is
complete. Native fixtures and compiler checks do not establish GUI usability or
experimental cloning success.

Validation passed: eight core design fixtures and all93 native web tests, including
six new library-state/persistence groups and two stable-role groups. The web WASM
compile passed. Cases cover wrong independent targets, equal-content/historical
source bindings, unselected unfinished preparations, per-construct raw edits,
pending-source persistence, role tampering and explicit version migration. No
browser checks or preview refresh were run for this checkpoint.
