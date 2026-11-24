# Original-source overlap assembly recipes

`core::assembly_recipe` compiles explicit preparation rows to existing typed
`DesignPlan` operations. It introduces no new scientific operation or plan schema.
The same compiler is shared by the browser worker and MCP host.

`Recipe` supplies an ordered set of `FragmentSpec` rows, explicit overlap lengths,
linear/circular output topology and `ProfileOptions`. Each row identifies a
stable original `PlanInput.id`, a distinct role ID, one preparation and an optional
reverse orientation. Repeated copies of the same original are allowed through
separate roles and remain separate ordered fragments.

Preparations are:

- `WholeLinear`: uses the entire actual original linear DNA document. A circle
  cannot become linear through this selection; no reference extraction is
  presented as prepared material.
- `ExactPcr {forward, reverse}`: emits the existing `Pcr` operation, preserving
  actual ordered primers, annealing domains and tails. The verifier freshly
  enumerates competing geometries and withholds ambiguous preparation.
- `Digest {enzymes, fragment_index, repair}`: emits the existing complete-cleavage
  `Digest` operation, selecting an explicit zero-based fragment in its existing
  top-reference cut ordering. Only enzymes supported by that operation are
  accepted. Optional `repair` adds a separate `RepairEnds` operation with explicit
  `FillFivePrimeOverhangs`, `RemoveThreePrimeOverhangs` or `Blunt` policy. `None`
  retains actual cohesive ends. Cleavage is the existing assumed-complete model,
  not a new methylation or condition model.

The compiler appends `ProfiledOverlap` with the declared ordered sources and
orientations. A reverse flag applies to the prepared fragment at assembly; it
never rewrites or reverses the original source document. `PreparedRow` records the
original ID/revision, preparation step, optional repair step and final prepared
source. Annotations and source evidence are carried by existing verifier
transformations, rather than a second coordinate-mapping implementation.

`compile(inputs, target, recipe)` validates structure and immutable source
revisions, then returns `Compiled {plan, rows, assembly_step}`. Compilation does
not assert that the plan works. The target is supplied independently and bound
by its actual revision; it is never generated from an assembly product or used
to choose PCR geometry, digest fragment, source order or orientation. Equal DNA
or even equal snapshots do not establish independence or invalidate a legitimate
explicit role choice. Hosts retain document/event identities and require an
explicit independent target selection; those identities are not inferred from
sequence content.

`preview(inputs, target, recipe)` runs the existing full verifier once. It returns
the compiled plan, complete verification evidence, and fresh prepared documents
in original orientation. If a preparation fails, `prepared` contains only the
available **leading prefix** of rows, never filtered successes with missing
preceding rows. An advisory overlap profile is calculated from fully available
prepared fragments using the existing profile engine. It is useful even when an
independent target differs. A profile failure is explicitly recorded in verifier
findings and returns no profile. Reference geometry does not override physical
substrate checks: staggered/nicked/blocked preparations may require additional
explicit supported operations before assembly. A profile can therefore be valid
while the physical plan remains incomplete. No repair is silently applied.

Fresh profile scans can occur twice: once in the verifier and once to return the
structured junction report. Both use the declared finite profile budget. Neither
published homology recommendations nor a consistent target comparison establish
reagent compatibility, kinetics, experimental success or exhaustive alternative
assembly uniqueness. Imported reference-only or physical-intermediate metadata
cannot grant new chemical authority.

`source_preflight(&[(&str, &Molecule)], target, recipe)` is a borrowed guard for
hosts before cloning originals. `preflight` accepts existing `PlanInput` objects.
`preview_revision` streams a bounded canonical hash over the full original input
bindings, actual target and recipe, including repair/profile choices. Readonly
`Preview::validate`/`validate_binding` check bounded models, prepared source
mapping, profile revisions/options and request context without replaying biology.
A decoded preview never restores private worker state or authorizes an edit;
candidate execution still requires fresh plan verification.

Bounds are 128 inputs/rows, at most 257 typed steps, 8 MiB total input-plus-target
DNA, 16 MiB source DNA counted over row reuse, and 32 MiB canonical request/output.
Existing borrowed native-model limits apply before hashing/cloning. Exact-PCR
annealing domains are 1–1,000 concrete bases and tails 0–1,000 bases. Aggregate
preparation scan work is charged before execution: PCR `4*N*(F+R)` and digest
`2*N*sum(recognition lengths)`, at most 100 million comparison units. Existing
PCR pairing/product bounds and digest/profile/verifier bounds remain in force.
Every junction has an explicit positive overlap up to 10,000 bases, with total
homology up to one million and no more than the declared profile comparison
budget. Original source attachments remain in the native plan, while host
history retains original document events.
