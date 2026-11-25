# Explicit review of observed read differences

`read_review` records human or caller decisions over one saved construct-owned read check. Preview and commit each freshly replay the exact original reference and ordered read events once under the current read-validation engine. Imported historical reports never supply current call eligibility. No original reference, chromatogram, FASTQ calls, quality evidence or document head is overwritten.

A request names the saved check and proposed reviewed reference, and supplies decisions with explicit reasons. Evidence keys identify a reference site or an insertion report index. Retaining the reference and leaving evidence unresolved are record-only decisions. They do not establish agreement, fill uncovered consensus calls, or create a new target.

A selected base includes the read index, current source-event call index (`ReadColumn.call_index`) and **expected concrete base**. The fresh report must contain that exact eligible call at the selected reference position, with the same expected A/C/G/T call. A selected deletion requires an actual eligible deletion observation from the nominated read at that site. A selected insertion includes its expected read, reference boundary and complete concrete sequence; all must match the fresh report and every eligible original insertion call. Changed placement/index meaning causes rejection instead of silently applying another allele. Deliberately choosing an eligible allele from conflicting observations is permitted with an explicit warning; the full conflicting report remains retained and is not relabeled consensus.

Duplicate evidence keys and contradictory insertion choices at the same boundary are refused. Circular boundaries zero and sequence length share one boundary for conflict detection. Changes use original reference coordinates in descending order. At a shared boundary, the original base is replaced/deleted before an insertion is made there. `Molecule::replace` preserves supported feature-coordinate handling, affected-feature review flags and source-bound citation behavior; unsupported coordinate-bearing evidence remains an explicit error. Empty resulting DNA is refused by the native model.

Changed targets are restrictive reference-only documents, carrying the saved-check hash, exact reference/read bindings, explicit edits and decision provenance. Physical descriptors are retained as ancestral evidence; these outputs are not intact duplex substrates or experimental confirmation. A commit derives a new document with the original reference **and all observed read events** as computational parents, including reads whose calls were not chosen. No target is created when decisions make no DNA change.

Every saved `read-review-v1` record contains schema version 1, the request/decisions, exact original reference/reads, one complete fresh report/options, actual edits, limitations and optional historical target event. Installation contexts bind the current complete workspace and request. Saved evidence identity is separately stable across retries. Repeating a commit still replays current science, then reuses an identical current named target/record. If a previously produced target has been edited, its head is preserved and a separate newly derived target is created; subsequent retries reuse that result. Earlier review records and original source history remain available.

`preview_revision`/`commit_revision` and borrowed `source_preflight` precede host copies. Limits are 256 decisions, 1024-byte reasons, a 256-byte name, 8 MiB source DNA, 100,000 replacement bases, 128 review attachments and 32 MiB per record. Existing native read/model/workspace caps apply. Conservative source-copy/annotation work is charged before repeated edits. Insertion calls use a bounded index, avoiding a per-inserted-base scan over all alignment columns.

`load`/`list` expose historical evidence. Readonly preview/commit validation checks context, bounds, explicit choices, report identities, restrictive target markers and actual retained source-parent events without rerunning alignment. Loading or decoding cannot authorize experimental claims. Hosts must compare request revisions before installing a returned workspace, and display caller decisions separately from the native scientific report.

Five independent fixtures cover a literal nominated allele, original-parent preservation and retry deduplication; retain/unresolved record-only decisions; explicitly selected conflicting evidence; full insertions/deletions and shared-boundary ordering; circular seam conflicts; changed expected alleles/positions; and edited-output preservation with a reusable new derivation. No read/alignment algorithm, consensus rule, scientific operation or verification engine version changes.

## Connected review controls

The construct's saved-check panel opens review in a central dialog. Fresh replay is
required before evidence choices become available. Sites and insertions page in
bounded groups; explicit Add/Update/Remove actions retain the request draft. Notes
are optional in the UI: blank notes generate factual choice-specific provenance,
while custom text is retained. Changing a choice refreshes an automatic note instead
of retaining a statement about a previous allele. Unselected sites remain unreviewed.

Preview is separately reconstructed and read-only. Saving uses a cancellable private
worker and one IndexedDB current/previous recovery transaction; no partial review or
new reference is installed on failure/cancellation. Record-only saves preserve the
pinned construct; a new derived reference becomes the explicitly checked construct
and requires fresh read validation. Recorded decisions reopen from the original
checked construct, seeding editable choices without restoring fresh evidence.

Five core fixture groups, two binary protocol cases and two real MCP SDK groups
passed. They include record-only decisions, literal allele changes, insertion/
deletion ordering, exact source parents, unsupported evidence rejection, repeated
commit deduplication and preservation of an edited prior output. The full128-tool
stdio smoke and final web/worker WASM compiler check also passed. Browser checks
remain deferred; source/native checks do not establish GUI acceptance.
