# Ordered original-source ligation recipes

`core::ligation_recipe` compiles declared preparation rows into existing typed
operations and previews them against an independently supplied target. No new
scientific plan operation or report layout is introduced. Existing overlap recipe
DTOs retain their schema31 layout.

`Recipe` contains a name, ordered `Fragment` rows, product `circular` topology and
explicit `allow_nicked` policy. Each fragment contains a shared
`assembly_recipe::FragmentSpec` plus optional terminal `EndTreatment`.
Preparations remain `WholeLinear`, actual `ExactPcr`, or selected complete
`Digest` with optional separately declared end repair. Shared borrowed preparation
guards and `compile_preparations` serve both overlap and ligation compilers;
there is no synthetic overlap recipe and no duplicated PCR/digest algorithm.

`compile(inputs, target, recipe)` preserves full original `PlanInput` molecules
and their revisions. It emits existing `Pcr`/`Digest`/`RepairEnds` operations,
explicit `TreatEnds` operations where requested, then a final ordered `Ligate`.
The shared preparation mapping retains source IDs, distinct roles, preparation
and repair step indices. A ligation `Row` additionally records its treatment step
and final prepared source. Reusing an original under different roles does not
deduplicate or reorder material. Reverse orientation is applied to the actual
prepared duplex at assembly, with strand/end transformations handled by the
existing chemistry engine.

`preview(inputs, target, recipe)` runs the verifier once, receiving an internal
fresh physical trace from that same reconstruction. It returns:

- Full compiled plan and verification findings/products/annotation evidence.
- `prepared`: the available leading row prefix of reference documents, before
  assembly orientation. Failed rows never cause later successes to shift indices.
- `physical`: fresh physical summaries aligned with that prefix, in assembly
  orientation. Missing state is explicitly `None`.
- `junctions`: adjacent upstream-right/downstream-left end evidence, including
  a separate last-to-first junction for circular closure. These are endpoint
  descriptions; they do not independently claim compatibility or sealed bonds.
- `product_physical`: the actual final duplex summary only when reconstructed.
  Its retained nick positions and closure state remain visible even if the
  separately bound target is wrong.

PCR/digest/repair/treatment physical states come directly from freshly tracked
verifier intermediates. No `duplex_model` or imported producer metadata is parsed
to restore authority. For untreated `WholeLinear`, the existing intact-reference
constructor is the only starting assumption: free ends are **unknown**, not
implicitly phosphorylated or assigned accessible 3′-OH groups. Single/mixed-strand
and imported intermediate markers restrict that constructor. Those preparations
can therefore remain incomplete. `Phosphorylate` changes 5′ chemistry only and
cannot resolve an unknown 3′ terminus. This recipe does not add a source-chemistry
declaration operation; use supported fresh PCR/digest preparation when known
terminal chemistry is required.

Ideal exact PCR produces 5′-OH/3′-OH blunt material. An explicit phosphorylation
step can make both strands sealable in this model. Omitting one side's phosphate
can leave a nick only when at least one bond forms and `allow_nicked` is enabled.
Neither side sealable is rejected even under nick permission. End repair and
terminal treatment are explicit operations; internal nicks are never silently
repaired. Digest geometry, cohesive polarity, complementary overhangs, strand
orientation and circular closure are checked by existing ligation chemistry.
A consistent plan remains conditional on those declared assumptions, not kinetic
mixture uniqueness, experimental yield or biological success.

The target is an independently supplied document with its own revision. It does
not choose source rows, digest fragments, PCR geometry or orientations. Equal
DNA is permitted; hosts pin source/target document events and explicit roles
rather than infer independence from sequence content.

`source_preflight` borrows original ID/model pairs before host copies.
`preview_revision` binds full originals, target and recipe by bounded streamed
canonical JSON. `Preview::validate`/`validate_binding` check shape, model bounds,
row/treatment source mapping, endpoint-summary correspondence and request context
without rerunning chemistry. Decoded evidence is readonly and cannot authorize
an edit or replace fresh verification.

Shared preparation limits apply: 128 rows/inputs, 8 MiB total input-plus-target
DNA, 16 MiB role-weighted source DNA, 100 million precharged preparation comparison
units, existing finite PCR-pair/digest/verifier bounds and a 32 MiB canonical
request/output ceiling. Ligation adds at most one treatment per row, so there are
at most 385 operations. Physical summaries/nicks and aggregate output are bounded.
No external software, network service or target-guided search is used.

The [literal three-piece example](../fixtures/ligation-demo.md) pairs the worker protocol fixture with explicit preparation and target instructions.
