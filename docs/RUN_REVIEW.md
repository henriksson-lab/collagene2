# Saved-run review preparation

`run_review::prepare(workspace, run_id)` resolves the retained run's separately declared intended target and original inputs by exact document, event and revision bindings. It preserves the complete saved operation DAG, options and intended orientations in a `transfer::Candidate`. Current document heads are neither substituted nor moved. An equal-DNA document or event is not an acceptable replacement.

`PreparedReview` includes the original target/event, ordered input bindings, source workspace fingerprint, run attachment hash, saved assessment, current verifier engine and explicit historical-engine flag. The recorded assessment is documentary only. Preparation does not verify, reconstruct physical material, select an alternative product, create a target, install a candidate, export a protocol or update a biological head. `validate_for_workspace` recomputes the full documentary receipt and rejects modified or stale preparations, including workspace changes unrelated to the selected run.

Current-format records can retain a documentary assessment produced by a different scientific engine. A plan containing saved-material reconstruction retains those exact origin recipes. Preparation does not translate material models: current ordinary verification can refuse origins whose scientific engine differs. The UI must show this limitation and require fresh verification before making current acceptance claims; a current candidate container is not current scientific authority.

The helper has a 32 MiB streamed serialized receipt bound in addition to existing native workspace, run, candidate and shared-origin limits. It never needs to serialize the entire workspace into a second buffer. Worker, web and MCP adapters call this same native helper.

Focused fixtures cover a real retained PCR with subsequently edited stock and intended-target heads; equal-snapshot event aliases and missing original events; candidate/stale-context tampering; and a structurally valid historical-engine archive whose status remains documentary. No GUI or downstream physical migration is claimed by this module.

## Saved-run UI

The run-details panel embeds a dedicated review section. Its explicit action sends `PrepareSavedRunReview` to a disposable worker, validates the receipt against the exact workspace, then sends the preserved plan and retained independent target to ordinary current `Verify`. Findings and reconstructed products remain read-only; retained ligation alternatives use the same current evidence viewer. Operation/product lists are paged. No candidate is installed or biological document moved.

The section reuses `AssemblyOutputs` for current HTML/PDF and oligo exports. Its capture callback checks workspace epoch, full workspace fingerprint, exact run attachment, prepared receipt and current report identity before work and again at native file acceptance. Export completion follows the existing acknowledged native file path. Exports are allowed for diagnostic invalid/incomplete reviews with their actual current status; they do not promote them to accepted material. Closing or invalidating the panel terminates its cached worker; nested output workers have their own cleanup and guards. Generic editor handoff remains deferred rather than silently changing the active target.
