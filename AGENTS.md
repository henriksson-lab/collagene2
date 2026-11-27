# Repository rules

- Fix or delete old, obsolete code. Never add compatibility scaffolding.
- Prefer temporarily broken code over scaffolding, including during implementation and refactoring. Update the implementation, callers, and tests directly; intermediate builds do not need to stay green. Complete the change and required checks before declaring it finished.
- Keep one current implementation and one current internal plan/protocol representation. Update callers directly when changing an interface.
- Do not preserve superseded development formats through legacy branches, adapters, aliases, fallback decoders, or migrations. Reject unsupported internal versions explicitly.
- Remove tests whose only purpose is supporting obsolete implementations or encodings. Preserve meaningful scientific correctness, data integrity, resource bounds, and current-format validation tests.
- Intentional import/export support for external file formats, including original Collagene files, remains product functionality. It does not justify compatibility layers around internal code.

## Verification

- Do not run browser, CDP, screenshot, or application HTTP checks until in-scope code is complete, as requested by the user.
- Run necessary compiler checks in the background while continuing implementation. Use the existing shared Cargo cache; do not introduce additional caches or automatic release-build verification.
