# Saved material digests

Saved digests use the same target-free material derivation records as reusable exact-PCR products. There is one processing path: select an exact native source event, declare a catalog-pinned digest and fragment index, reconstruct it natively, and retain the resulting Derived material. This processing does not invent an independent fragment target.

For an existing saved run, construct `material_derivation::Request` with `SourceOrigin::Saved(MaterialOrigin::Run(capsule))` and `Derivation::CatalogDigest { selection, options, fragment_index }`. An already derived source uses `MaterialOrigin::Derived`. The request also declares the output name and whether to retain newly declared primer orders. Installation retains the canonical native record, product event and source dependencies together.

The browser digest workflow uses `InstallMaterialDerivation`. Its source picker resolves the exact installed event to a Derived capsule. `Preparation::FromDerived` emits an input-bound `ReconstructMaterial` step before the selected downstream preparation. Original source recipes, primers, cut choice and event pins remain intact; a saved fragment is not flattened to an ordinary sequence input.

Fresh replay retains actual top and bottom strands, strand offsets, free-end chemistry, nicks, positional modifications and inherited incompleteness. A reference contour or serialized duplex summary does not supply physical authority. PCR on staggered material remains restricted to supported geometries within its known reference; unsupported cases refuse explicitly. Downstream cloning retains its independently declared target and intended path.

The Run and Derived families share one iterative ancestry ledger. It charges wrappers, source models and opaque candidate/record bytes before decoding or cloning, under aggregate depth, byte, base and step bounds. Native ownership checks bind retained event/document/revision identities. Equal DNA, a renamed document or a producer-label string cannot replace the actual source event.

Oligo orders include only explicitly declared new preparations. Historical primers remain provenance and are not ordered again. HTML/PDF protocols describe the retained stock and new steps; ligation output preserves intended and bounded alternative paths.

Public and private selected-document publication follow the generic material-record index and the same bounded retention queues as native runs. The closure retains exact source events, owning runs, nested derived records and required evidence. Unrelated private documents, later private heads and editable drafts are excluded. Canonical records and source geometry remain freshly replayable after archive roundtrip.

See [material derivations](MATERIAL_DERIVATIONS.md), [controls](MATERIAL_DERIVATION_CONTROLS.md), [publication](MATERIAL_DERIVATION_SHARING.md) and [MCP](MCP_MATERIAL_DERIVATIONS.md). Current fixtures cover literal sticky/blunt geometry, modification ancestry, downstream PCR and assembly, incomplete origins, wrong-event refusal and publication replay.
