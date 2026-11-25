# Manual cloning plans and candidate edits

The Cloning plan panel now offers **New manual plan** and **Edit steps and sources**.
This reuses the existing native operations, independent target, candidate storage,
verification, product/junction review, retained history and HTML/PDF outputs.
A new draft can start empty; it is saved as a candidate only after explicit Apply
with at least one structurally valid step and exact source bindings.

## A connected everyday workflow

A [loadable synthetic example](../fixtures/manual-adaptor-demo.md) provides all
four GenBank records and exact field values. To clone an annealed adaptor into a backbone:

1. Select a separately supplied intended target and open a new manual plan, or edit
   the candidate from the existing oligo-annealing panel.
2. Add the two original oligo inputs and original backbone, choosing exact document
   events. Insert oligo annealing with explicit strand order, pairing threshold,
   selected candidate and terminal chemistry.
3. Insert a backbone digest. Add ordered ligation using the digest and annealed
   products as earlier-step sources; choose orientation, topology and nick policy.
4. Apply the declarations, then Verify. Review actual products and junctions before
   finalizing named outputs or exporting protocols. Failed checking remains visible;
   editing supplies no positive scientific evidence.

PCR fields include names, annealing domains and separate 5′ tails. New PCR fields
start empty and require entered domains. Digest fields expose enzymes and a
zero-based fragment index. End treatment, repair, single-base tailing, annealing
and ordered ligation have typed controls. Existing complex method operations keep
their preparation/settings and allow source rebinding; their detailed preparation
controls remain in the existing method panels. Bound analytical selections can
need fresh inspection after edits.

## Dependency edits and exact identities

The core `plan_editing` helper validates named input references and strictly earlier
step dependencies. Moving a step remaps producer indexes so consumers still refer
to the same producer. Removing a referenced producer is refused; remove or rebind
its consumers first. Insertion shifts existing dependency indexes, while the new
operation selects pre-insertion earlier producers. A rejected edit leaves the
original step vector intact. These helpers evaluate no DNA or chemical state.

The editor holds input models through shared references while field edits affect
operation drafts. Long lists page ten steps/inputs and six source references.
Replacing a step and applying the whole draft are separate explicit actions.
Cancelling discards the local draft. Existing saved candidates remain available.

The host preserves exact preferred source events. Other original inputs require
an explicit event choice, including historical events with the required revision;
content equality never selects their document ancestry. Current DNA heads remain
available and additional historical choices have an explicit bounded-list notice. Added inputs use a chosen document/event/revision and a unique input ID.
Original input revisions and the independent target stay fixed. A source is
resolved against its retained document history before Apply and fresh Verify; a
newer current head never replaces the chosen event. Stale workspace,
active-candidate or target changes reject the commit. The new candidate gets new
source preferences and clears previous verification/product-finalization state.
Fresh native reconstruction alone establishes current scientific material.
Candidate event preferences are session state. Reopening a portable candidate can
require fresh source choices; retained runs preserve exact events in their native
history. Unapplied editor field changes are local and Cancel discards them.

## Native validation and compatibility

The independent adaptor fixture supplies top `AATTGCTACG`, bottom `AATTCGTAGC`
and circular backbone `GGAATTCTTACG`. A six-base overlap produces two ten-base
strands with AATT ends; actual EcoRI digestion cuts backbone top/bottom at 2/6.
Ordered ligation yields the frozen 22 bp target `AATTCTTACGGGAATTGCTACG`, accounting
for 44 source nucleotide calls. Valid step movement preserves the physical result;
invalid ordering, referenced deletion, wrong enzyme and wrong independent target
cannot inherit a successful result. Candidate/archive/replay retains all original
events and excludes the target from physical product parents.

Graph edit fixtures independently track producer identities and preserve primers,
tails, orientations and operation settings through valid edits and atomic failures.
Scientific checking remains the existing deterministic verifier. Fifteen focused
native cases pass, as does the full application/worker WASM build. Port8082 serves
the new immutable snapshot with compatibility workers; no application browser or
HTTP smoke checks were performed.

No operation/report/wire fields or runtime dependencies are added. Worker67,
candidate15, plan13, science18 and MCP165 remain the current interfaces. Other
unfinished method chemistry, richer preparation controls, durable graph-only
review, bulk/platform delivery and final GUI acceptance remain active. Browser
checks stay deferred; this increment makes no rendered-layout acceptance claim.
