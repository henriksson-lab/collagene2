# Source-library edits and recovery

The existing construct library supports reviewed bin insertion, removal and
reordering. These operations change editable design intent across the library;
they do not edit source DNA, rewind document heads or establish a cloning result.

## Explicit transformations

Reordering uses stable bin IDs and requires the complete permutation. Removing
the last bin and unchanged permutations are rejected. New bins require explicit
source alternatives and an explicit source/preparation choice for every existing
construct. Fresh captures use exact retained document events and original source
coordinates. They begin as whole-source preparations; choosing a contribution
interval does not invent a PCR product or physical extraction.

Users can copy an existing immutable source/preparation choice or capture an open
DNA document. Pool constructs also require explicit material choices. The add
form pages sources and constructs and offers explicit apply-to-all actions;
there is no automatic first-alternative assignment.

Surviving bin, choice and construct IDs remain unchanged. New identities consume
the monotonic library counter. Preparations, named primers and literal tails move
with their source bins. Directed bin adjacency identifies junctions: surviving
edges retain their declarations, while changed edges require review. New overlap
edges leave the profile uncomputed. All previous raw text and removed settings
remain in recovery. Pool candidate declarations are cleared after topology edits.

## Review and atomic publication

An opaque in-process proposal retains exact before/after drafts, a review revision
and per-construct impacts. The preview identifies primer neighbors needing review,
removed preparations, retained/removed/new junctions, overlap changes, source-role
order and candidate changes. Restorations additionally identify disappearing and
returning constructs. All impact rows are paged.

Applying requires the unchanged complete draft, workspace epoch and native
workspace fingerprint, plus acknowledgement of that exact proposal. The host
validates both full draft models and every captured source, including inactive,
removed and working choices. Equal DNA in another event is not a replacement.

The new active draft and immutable recovery record are staged in one private
workspace. A single cancellable IndexedDB transaction installs that frame with
the previous recovery workspace. No earlier autosave can install a lossy
intermediate. The same retained Rust worker encodes both frames; cancellation
closes pending encodes or aborts the transaction. Token-owned commit leases prevent
old cancelled completions from unlocking a subsequent write. Root-owned application callbacks outlive the dock or modal that initiated the
write. A durable successful commit reaches the application even if that panel
closes immediately afterward.

## Exact raw recovery and restoration

Recovery uses bounded raw serialization, preserving the distinction between
working fields and saved construct settings. The ordinary normalized encoder is
not used for these snapshots. Application draft autosave uses the same raw encoder,
so restored unsaved text is not silently normalized by a subsequent save.

The native `assembly-library-history-v3` record retains exact before bytes,
persisted-active/raw/after hashes, proposal and source-workspace revisions, source
event pins, typed impacts and allocation counters. Its Postcard attachment avoids
expanding raw JSON into integer arrays. Native validation checks structure,
hashes, retained event membership and allocation bounds; complete WebDraft
validation remains a separate host responsibility. Imported records are historical
editable data and confer no scientific authority.

History supports bounded listing, metadata/impact/source review and exact JSON
download. Restoration first validates the full recovery draft and every cached
source against the current workspace. It is another reviewed edit: current raw
work is saved before replacement. The same library identity is required, and the
allocation counter becomes the maximum of current, recovered and retained history
values. Previously allocated IDs cannot be reused. Old candidate choices remain
declarations requiring fresh native verification, never restored results.

## Bounds and MCP

Existing limits remain12 bins,64 choices per bin,256 total choices and256 selected
constructs, with32 MiB raw draft/request bounds. Restore impact unions permit512
construct records and the before/after source union permits536 exact event pins.
Recovery retains at most128 records and32 MiB aggregate history; a capacity failure
rejects the complete edit without truncation. Native archive and attachment limits
also apply. History pruning/retention management and broader cross-library recovery
remain delivery work.

MCP adds `list_assembly_library_history` and `open_assembly_library_history` over
immutable workspace references. Open returns metadata and exact recovery JSON as
separate artifacts. These149 registered operations include no arbitrary filesystem
path, remote message, scheduling or inventory capability. No new scientific
operation, worker selector or native schema is introduced: worker59, candidate12,
plan10, science15 and draft3 remain unchanged.

## Validation

Passed135 focused native checks: seven recovery-record fixtures, twelve native
workspace fixtures, all114 web state/model fixtures and two real MCP SDK groups.
The final WASM check passed for application and worker. New cases cover exact
unsaved working/saved-row differences, inactive source membership, stale draft/
frame/ack rejection, stable identity and explicit all-row additions, restore impact
unions, counter highwater/exhaustion, and paired metadata/raw artifact fidelity.
Active working pool changes require their own explicit material selection.

Browser, visual and performance acceptance remains deferred until the in-scope
implementation is complete. Port8082 retains the previous checked reaction-matrix
snapshot; no new browser checks or full preview rebuild were run here.
