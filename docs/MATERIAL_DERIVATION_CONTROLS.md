# Reusable PCR and digest derivation controls

`web/src/material_derivation_model.rs` contains pure model code registered for native fixtures and browser controls.
It prepares explicit native `material_derivation::Request` values. It performs no
browser transaction, file write, visible document import or successful-install
notification. The PCR and digest dialogs use the existing cancellable atomic recovery transaction before visible apply.

`SourcePin::capture(event, epoch, source, origin)` retains the exact document,
history event and revision, workspace epoch, and optional closed `MaterialOrigin`.
An equal sequence under another document/event is rejected. Changing or removing
an origin invalidates the pin; a saved origin never silently becomes an original
reference. Native source/record checks remain authoritative.

`PcrChoice::capture(pin, forward, reverse, index, Arc<PcrDocument>)` retains both
complete named oligos, selected geometry including its orientation and primer
identities, and the exact analytical preview revision. `validate_current` checks
all those fields before `request(workspace, name)` builds an explicit
`ExactPcr` derivation. No independent target is supplied or invented. A preview
remains reference evidence until the separately installed native derivation can
be reconstructed from its retained record. The native initial profile still
requires one actual geometry; choosing a display row does not override that rule.

`DigestChoice::from_preview(pin, molecule, catalog, options, preview, index)`
validates the generic native digest preview context and determined-fragment choice.
The lower-level `capture(pin, catalog, options, inspection_revision, fragment)`
retains the exact catalog pin, complete scan limits, inspection context and
selected two-strand geometry. The caller must first validate its native inspection
response. Current-state checks invalidate edited catalog/options/fragment choices.
`request` builds a `CatalogDigest` derivation; actual cleavage is recomputed.
The origin type also admits derived PCR/digest documents for repeated operations;
there is no Run-only cast or reference fallback in this model.

Both choices return `PendingInstallation`, with getters `request()` and
`context_revision()`. `validate_workspace(workspace, epoch)` pins the complete
original workspace fingerprint. `validate_received(installed, decoded)` checks
the native context/source/archive receipt and equality with the actual received
archive. This is **receipt validation, not transaction completion**.

The integrated UI rechecks the current source/choice and workspace before worker
dispatch, before persistence, and before visible apply; it then atomically stores
current and previous recovery frames using existing cancellation/rollback rules.
Only after that transaction succeeds does the integrated panel open the product and report success.
Cancellation or failed persistence must leave the reference preview and original
workspace unchanged. Findings and inherited incompleteness remain visible.

Five native fixture groups cover source/epoch/oligo/geometry invalidation,
validated receipt versus uncommitted workspace, explicit derived origins, and
catalog/limits/fragment invalidation, and actual EcoRI-opened circle followed by a
Derived-origin SmaI digest using the generic native inspection. Source is formatted; parent native/WASM checks validate this integrated phase.

The ordinary PCR output checkbox now explicitly opens a reusable product, rather
than promoting the analytical reference preview. `Request.save_primers` binds the
optional current primer collection to the same native installation/archive;
primer-only output retains its separate existing route. Generic digest inspection
accepts Run, Digest and Derived origins and sequentially opens fresh products
through that same installer. Assembly uses `FromDerived`; draft13 is required for
that origin and unknown future14 is rejected. Advanced run-only tools continue to
refuse unsupported origins explicitly.
