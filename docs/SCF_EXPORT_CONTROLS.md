# Reviewed SCF export controls

The chromatogram's compact **Export trace…** menu offers regeneration of the displayed actual trace as SCF3.00. The existing original-attachment download remains separate. There are no additional permanent export toolbar buttons.

Choose smallest lossless, explicit8-bit or explicit16-bit precision, then **Prepare export review**. The binary analysis worker reconstructs and validates the unchanged trace attachment and explicit orientation. It returns prepared bytes plus the documentary-loss inventory; this step does not write a file. Negative/unrepresentable signal, fabricated/edited trace fields, missing quality vectors and ABIF scalar scores are refused by the strict native profile rather than converted into invented trace evidence.

The review shows actual file size, precision, orientation, original source hash, warnings and omission categories/counts. Any omissions require a separate checkbox acknowledgment before **Save reviewed SCF…**. This acknowledgment permits only the displayed documentary omissions; signal/calls/peak/quality loss remains unsupported.

Prepared results bind the host workspace epoch and exact document/event identity, full current molecule revision, source-file hash, plotted orientation, selected precision and worker context. A document, orientation or precision change closes active preparation and clears the review/acknowledgment. The worker reply must match request ID/context, source/orientation and precision, and pass the native report's documentary validator. Closing/unmounting the controls invalidates late replies.

Save rechecks the same scope/acknowledgment before reserving export and again at native export admission. The panel reports success only after the existing acknowledged file-export helper completes; native receipts show the real managed export basename centrally. It does not fall back to a browser Blob on native failure or claim success for queued work. A stale panel cannot display a late completion against another document.

`scf_export_model` provides native-reachable pure guard fixtures for explicit omission acknowledgment, changed document with identical attachment, separate equal-content document/event/workspace identities, and orientation/precision invalidation. Register the model with `cfg(any(target_arch="wasm32",test))`; the view is WASM-only. Core binary roundtrip/layout fixtures remain in `crates/core/tests/analysis/reads/scf_export.rs`. These controls do not claim vendor-application compatibility or biological verification.
