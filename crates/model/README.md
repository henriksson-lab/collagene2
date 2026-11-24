# collagene-model

Shared MIT molecular document primitives extracted without behavior changes from Collagene's core: `Alphabet`, `Topology`, validated `Sequence`, `Location`, `Feature`, `Molecule`, `Error`, `complement`, and the existing size bound. Core continues to reexport these exact shared types through `core::model` and its existing top-level names.

There is no network, filesystem, UI, scientific algorithm dependency or parallel duplicate molecule model here. Serde normalization/validation, location extraction/edit behavior, enum layouts, field order and `Molecule::revision()` remain unchanged. The existing `json-schema` feature enables optional Schemars 1.2; the core forwards its feature to this crate.

Literal pre-extraction JSON, raw Postcard and revision fixtures live in `fixtures/`; tests also cover independent biological expectations, invalid calls/locations, recursive decode bounds, document constraints and retained edit-review semantics. Model extraction introduces no document, worker, plan, workspace or scientific schema/version change.

`bounds` provides the shared borrowed annotation/tree/string resource guards before cloning, recursion, hashing or serialization. The existing recombination guard paths reexport these same functions. Their aggregate limits remain unchanged; they do not grant molecular validity or alter serialized models.
