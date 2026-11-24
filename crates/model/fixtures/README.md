# Pre-extraction compatibility fixtures

Captured on 2026-11-24 **before** `crates/core/src/model.rs` became a reexport. A temporary native Rust harness called the existing `collagene_core::model` API, serialized the constructed values with the existing locked `serde_json` and `postcard`, and recorded `Molecule::revision()`. The original source's SHA-256 is retained in `original-model.sha256`. The new crate library is a byte-for-byte copy of that source.

These synthetic fixtures are deliberately small and independently specified: simple DNA with whitespace/case normalization; a circular molecule with a complement of an origin-crossing join, an optional qualifier, metadata inserted in reverse lexical order, and uninterpreted original evidence; and an ordered fuzzy/remote location tree. No real laboratory document or another application was used. The expected bytes and revision strings are fixed data; tests do not regenerate their expectations from the new implementation.

The Postcard files are raw model payloads, **not** new native/protocol envelopes. Existing envelope versions, checksums and schemas remain owned by the core crate. JSON files contain compact UTF-8 without an added newline.
