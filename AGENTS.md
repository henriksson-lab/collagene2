# Repository rules

- Fix or delete old, obsolete code. Never add compatibility scaffolding.
- Prefer temporarily broken code over scaffolding, including during implementation and refactoring. Update the implementation, callers, and tests directly; intermediate builds do not need to stay green. Complete the change and required checks before declaring it finished.
- Keep one current implementation and one current internal plan/protocol representation. Update callers directly when changing an interface.
- Do not preserve superseded development formats through legacy branches, adapters, aliases, fallback decoders, or migrations. Reject unsupported internal versions explicitly.
- Remove tests whose only purpose is supporting obsolete implementations or encodings. Preserve meaningful scientific correctness, data integrity, resource bounds, and current-format validation tests.
- Intentional import/export support for external file formats, including original Collagene files, remains product functionality. It does not justify compatibility layers around internal code.

## Verification

- Implement in coherent batches across the scientific core, schemas, transport, UI, callers, and fixtures before compiling. Prefer temporarily broken intermediate code to small changes made only to keep each step building.
- Reduce compiler and test invocations: inspect callers and enum matches statically, finish the batch, then run its focused native and WASM checks. Collect concrete failures and fix them together before rechecking. Do not rerun unchanged passing cases without a specific reason.
- Parallel workers contribute to the same integration batch; their individual patches do not each trigger a build. Keep Cargo package/feature selections consistent so rechecks reuse the existing cache.
- Do not run browser, CDP, screenshot, or application HTTP checks until in-scope code is complete, as requested by the user.
- Run necessary compiler checks in the background while continuing implementation. Use the existing shared Cargo cache; do not introduce additional caches or automatic release-build verification.
