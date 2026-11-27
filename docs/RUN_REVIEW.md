# Saved-run review preparation

`run_review::prepare(workspace, run_id)` resolves the retained run's separately declared intended target and original inputs by exact document, event and revision bindings. It preserves the complete saved operation DAG, options and intended orientations in a `transfer::Candidate`. Current document heads are neither substituted nor moved. An equal-DNA document or event is not an acceptable replacement.

`PreparedReview` includes the original target/event, ordered input bindings, source workspace fingerprint, run attachment hash, saved assessment, current verifier engine and explicit historical-engine flag. The recorded assessment is documentary only. Preparation does not verify, reconstruct physical material, select an alternative product, create a target, install a candidate, export a protocol or update a biological head. `validate_for_workspace` recomputes the full documentary receipt and rejects modified or stale preparations, including workspace changes unrelated to the selected run.

Current-format records can retain a documentary assessment produced by a different scientific engine. A plan containing saved-material reconstruction retains those exact origin recipes. Preparation does not translate material models: current ordinary verification can refuse origins whose scientific engine differs. The UI must show this limitation and require fresh verification before making current acceptance claims; a current candidate container is not current scientific authority.

The helper has a 32 MiB streamed serialized receipt bound in addition to existing native workspace, run, candidate and shared-origin limits. It never needs to serialize the entire workspace into a second buffer. Worker, web and MCP adapters call this same native helper.

Focused fixtures cover a real retained PCR with subsequently edited stock and intended-target heads; equal-snapshot event aliases and missing original events; candidate/stale-context tampering; and a structurally valid historical-engine archive whose status remains documentary. No GUI or downstream physical migration is claimed by this module.

## Saved-run UI

The run-details panel embeds a dedicated review section. Its explicit action sends `PrepareSavedRunReview` to a disposable worker, validates the receipt against the exact workspace, then sends the preserved plan and retained independent target to ordinary current `Verify`. Findings and reconstructed products remain read-only; retained ligation alternatives use the same current evidence viewer. Operation/product lists are paged. No candidate is installed or biological document moved.

The section reuses `AssemblyOutputs` for current HTML/PDF and oligo exports. Its capture callback checks workspace epoch, full workspace fingerprint, exact run attachment, prepared receipt and current report identity before work and again at native file acceptance. Export completion follows the existing acknowledged native file path. Exports are allowed for diagnostic invalid/incomplete reviews with their actual current status; they do not promote them to accepted material. Closing or invalidating the panel terminates its cached worker; nested output workers have their own cleanup and guards. An explicit Edit retained plan action hands a separate candidate to the existing DAG editor. It revalidates the full native receipt and workspace epoch, retains exact source and target events, clears cached reports and leaves viewer heads unchanged. Edits require fresh verification and create a separate run when recorded. Editing preferences are session-only; after restart, prepare the saved run again.

## Native editing acceptance

The native editing fixture prepares a real saved PCR after both its stock and
intended-target document heads have advanced. It keeps the review receipt
immutable, edits a separate candidate through `plan_editing::inserted`, adds an
explicit phosphorylation step, verifies against the retained original intended
target, and generates current HTML/PDF. Recording the new run uses the receipt's
original source and target events. The existing heads, history prefix, attachments
and original run assessment remain unchanged. After native archive restart, the
new run replays successfully and yields actual phosphorylated PCR material.

At admission the caller rechecks the original receipt against the current native
workspace. Equal-snapshot aliases cannot replace its exact source or target event.
An unrelated workspace edit invalidates it, as does a changed run attachment byte
pin even when the re-encoded record still has identical biological declarations
and assessment. Generic `record_run` separately accepts explicitly declared new
source bindings; the saved-review caller must retain these receipt checks rather
than treating equal DNA as permission to replace a historical event. The web handoff uses the same immutable review and exact event bindings. It does
not adopt a saved status or mutate the active working document.
