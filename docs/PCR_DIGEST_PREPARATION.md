# PCR followed by catalog-pinned digestion

`assembly_recipe::Preparation::PcrCatalogDigest` is an additional shared preparation for overlap and ligation recipes. Its fields are `forward`, `reverse`, `template_strand`, `selection`, `fragment_index`, and optional `repair`. Exact saved origins use the existing `FromRun` or `FromDerived` preparation wrapper and retain native source ownership.

The compiler emits existing operations against the original input or an explicitly reconstructed native source:

1. `Pcr { template: Input(original) or Step(reconstruction), forward, reverse, template_strand }`.
2. `CatalogDigest { template: Step(pcr), selection, fragment_index }`.
3. An optional, separate `RepairEnds` operation.

`PreparedRow.preparation_step` remains the first step in the preparation chain, here PCR. `prepared_source` identifies the final digest or repair step. A later end treatment in a ligation row remains another explicit operation. All intermediate products remain in the normal verifier/run history; no copied source interval is substituted for physical preparation.

The existing exact-PCR verifier retains alternative products and supports only its existing concrete-primer geometry. The existing catalog digest verifies selected enzyme identity and geometry against the embedded catalog, retaining actual top/bottom strand lengths, staggered ends and modification uncertainty. Complete cleavage remains a declared assumption. Added recognition sites can be supplied by explicit 5′ primer tails; selecting this preparation does not design those tails or assert that their cleavage produces the intended fragment. Internal sites can create additional fragments. The explicit fragment index chooses an actual determined fragment, never a target-guided matching interval.

Borrowed preflight validates primer lengths/alphabet, catalog selection and source models before compilation copies or revision hashing. The existing aggregate preparation-work cap is charged with both PCR binding work `4 × N × (F + R)` and a conservative digest scan `A × 2Σ(recognition lengths)`, where `A = N + forward_tail + reverse_tail` bounds one-turn PCR product length. The aggregate role-base budget includes both tails for this preparation. Existing engine-level enumeration, work and serialization limits still apply. This preparation does not change ordinary primer-order export's separate 200-base profile.

Four focused fixtures independently specify PCR-added BsaI sites, actual retained top/bottom geometry and overhangs, original input immutability, complete intermediate retention, optional repair chaining, overlap receipt mapping, wrong targets, additional internal sites, invalid fragment indices, tampered catalog pins, stale input revisions and aggregate work rejection. No new digestion or PCR algorithm, operation tag, scientific claim or external executable is introduced.

Standalone exact-PCR digest inspection accepts original references or generic
`MaterialOrigin::Run/Derived` with its exact owning workspace. The saved path uses
`preflight_from_material`, `inspection_revision_from_material` and
`inspect_from_material`, returning `MaterialPreview`. Worker108 transports the
same source/origin/request/workspace; MCP and the reviewed fragment picker use
that same contract. No old Run-only inspector remains.

Native Bottom frames retain real stored-strand calls beyond the displayed Top
reference. The independent13-Top/17-Bottom fixture retains17 template calls through
a target-free derivation and archive restart, then adds literal EcoRI tails and
checks the resulting physical digest ends. Inherited incompleteness and bounded
findings remain explicit. Inspection does not install material; changing primers,
catalog, exact origin event or source ownership invalidates fragment application.
Saved Derived mismatch-PCR and standalone OE digest inspection retain their
explicit unsupported diagnostics where their separate adapters are not connected.
