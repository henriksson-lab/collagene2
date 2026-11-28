# Controller and scientific owner integration

The current batch moves complete behavior into existing owners and updates callers
directly. App decreases from 2,869 to 2,413 lines. ProtocolController owns the
complete plan-verification lifecycle; RestrictionController owns saved enzyme
sets and digest actions; BiologyController owns antibody/RNA actions; and
ComparisonController owns comparison persistence and admission. No universal
callback container, forwarding App API or additional serialized state is added.

PlanState now owns its input-event bindings alongside the report, plan, target and
generation. Review found a preexisting concurrency error: a new verification read
the generation without advancing it, so distinct requests against unchanged
source events could publish in completion order. Each dispatch now invalidates
the complete evidence owner and captures the new generation. Completion first
rejects a superseded/disposed owner, then retains the existing workspace epoch,
exact target event/revision and complete source-head checks. Source changes and
restore also clear the complete evidence, including input-event bindings.

Preparation owns catalog lookup and operation counting. TrackedDuplex owns
summary, checked-digest preflight and private strand/index access. Their types and
implementations remain in the same existing files. Core, web, MCP and fixtures
use the methods directly; old free-function entry points are deleted. Multioperand
ligation, homology assembly and other stateless algorithms remain free functions.
Renderer State owns its five gesture/position/feature/draw operations in the same
file as its state and GPU lifetime; geometry and tessellation remain free.

Static comparison preserves 23 App initializer expressions, all six moved core
method bodies and all five renderer bodies after documented receiver/formatting
normalization. The verification-generation and complete invalidation fixes are
explicit exceptions to preservation. The
[accounting record](controller-owner-integration.json) distinguishes source
comparison from compiler/test evidence; token preservation is not runtime proof.

Native/WASM compilation and **107 focused cases pass**. All **2,509 baseline
case names / 33 harnesses** remain discoverable, with no additions, removals or
renames. The preceding `82937c1` checkpoint passed scientific and desktop CI;
exact-source CI remains required for this batch. Integration uses the existing shared Cargo cache and consistent package/features. Execute
changed fixture modules and affected native ligation unit case; compare
the complete executable discovery against the 2,509-case / 33-harness baseline.
Unchanged pure web-model cases cannot prove WASM-only controller callbacks and
are not rerun as a substitute. Browser, application HTTP, visual and installed
platform acceptance remain deferred. Scientific algorithms, hashes, supported
bounds, binary DTOs, schemas and engine stamps are unchanged.
