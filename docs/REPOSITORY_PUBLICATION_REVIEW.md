# Reviewed publication from a saved local catalog revision

The native remote UI must use the saved local library document ID and immutable library revision, not upload a working-editor archive. Existing `queue_push` and `queue_push_pinned` remain strict: an unsupported original SBOL graph is refused rather than silently shared or rewritten.

The `publication_review` API adds an explicit, local-process acknowledgement path:

1. `Repository::push_review_choices(id, expected_revision, remote)` reloads canonical state and inspects only structurally reachable sequence alternatives in that saved archive. Required sequences come from actual retained native claims. No first alternative is selected.
2. `prepare_push_review(id, expected_revision, remote, choices)` derives a fresh shared-library documentary subset review. `PushReview` exposes readonly selected references, original/subset hashes and omission counts through the shared review; no unrelated omitted identifiers, sequences or values are exposed. Preparation changes no outbox, local archive or network state.
3. The trusted host retains this opaque object. Only an explicit user acknowledgement calls `PushReview::approve()`. Neither `PushReview` nor `ApprovedPushReview` implements deserialization or has public mutable fields. A submitted review hash or reconstructed JSON receipt cannot manufacture approval.
4. `queue_reviewed_push(&approved)` reloads and checks the captured repository checkpoint, exact local catalog document and revision, original archive hash and native selection, configured remote, actual remote binding, and initial publication destination. It re-derives the approved projection from the immutable saved original and validates ordinary public document closure before encoding and enqueueing it.

An unrelated local placement of an already bound/mounted document is not interpreted as a remote move. First publication still requires an explicit destination folder binding. An intervening native library write, binding acknowledgement or repository change invalidates the review; prepare and review again. Approval cannot overwrite newer saved local work.

After enqueueing, the existing reliable outbox owns the exact public archive bytes, mutation ID, expected lab, native pin and remote base revision. Lost responses and restarts use ordinary receipt lookup and exact byte replay; no approval token or transient UI state is needed to resume. An explicit documentary-only update may save a new remote revision when its approved archive differs even if the native selection is unchanged. Identical bytes remain `NoChanges`; the strict push path retains its historical behavior.

The original full private workspace and SBOL bytes are unchanged locally. Subsets and counts-only omission receipts convey documentary scope, not verified SBOL derivation from absent originals or scientific/material authority. Native events, plans, exact source identities and physical capsule ancestry are preserved, and ordinary fresh replay remains necessary.

Focused native fixtures exercise an actual SBOL projection followed by saved PCR, real SQLite local storage and native authenticated dispatch, restart and a lost post-commit response, independent omission checks and fresh remote material replay, stale saved-head/metadata reviews, and explicit alternative-sequence publication with the same native head.


## Native controls and wire lifetime

The Databases accordion now exposes structured reachable-sequence choices and
explicit subset/hash/count review. It queues a saved publication without network
access; Pending pushes performs the separate deliberate resume. Prepared replies
carry the complete stored remote configuration for comparison with the selected
endpoint and lab. Dirty editor state, saved head, selection, configuration and
native checkpoint changes invalidate acceptance.

DSQ/DSR2 carries bounded metadata only. Legacy prefix commands still use version1
frames and receive compatible replies; new review commands/results require2.
The host keeps at most four opaque reviews behind random process-scoped IDs,
expires the oldest when replaced, and removes them on successful queue, discard
or credential-session reset. Cancellation is checked immediately before retaining
a token, so it cannot create an undisclosed retained review. No serialized review
hash or UI flag creates native publication authority. Lost durable queue replies
remain recoverable through Pending pushes after restart; approval IDs themselves
are intentionally not restart-persistent.
