# Ordered multi-fragment ligation

The Assemble / Cloning methods workflow includes ordered multi-fragment ligation alongside overlap assembly and TA/GC tailing. The active document remains the independently selected target. Choose original source documents in assembly order, prepare each with whole linear DNA, restriction digest (optional explicit end repair), or exact PCR, and choose its assembly orientation. Saved primer annotations may be explicitly copied into editable PCR name/domain/tail fields.

Each row separately declares no terminal treatment, 5′ phosphorylation or 5′ dephosphorylation. None preserves freshly reconstructed chemistry. Whole linear DNA sequence has unknown terminal chemistry; its origin/extension is not an assumed cleavage or experimentally established duplex end. Phosphorylation does not resolve unknown 3′ termini, repair geometry or infer host repair. Circular source documents require an explicit preparation.

Compute method evidence sends original source models and the full ordered recipe to a private worker. The shared core compiler/verifier reconstructs preparation, optional treatment, orientation and ligation. The UI shows actual oriented end overhang words and terminal groups, each candidate junction including circle closure, retained product nicks and the independent-target assessment. Missing physical summaries remain explicitly unavailable. Allow retained nicks is an explicit model option, never an implicit rescue. No unique mixture outcome or real reaction completion is claimed.

Review method candidate recompiles the original recipe into existing typed preparation, treatment and ligation operations. Prepared reference previews are not promoted into plan inputs or copied as physical authority. Source/config/order/target/epoch changes invalidate evidence; cancellation closes the current worker. Full method evidence JSON retains compiled source bindings, preparation and fresh physical/target findings. Existing candidate verification, protocols and native lineage remain the route for reviewing and retaining the plan.

The shared row form supports at most12 fragments. Initial Assemble mode is overlap; returning to the workflow preserves the user's selected mode and fields. This source-only integration has not run compiler or browser checks; those remain coordinated by the parent task. Mobile remains deferred.


The [three-PCR-fragment example](../fixtures/ligation-demo.md) independently
specifies an 18-base target, phosphorylation steps and a deliberate retained-nick
case. Native candidate capture preserves exact source document/event bindings,
including when different open documents have equal sequence content. Returning
to Assemble retains the chosen workflow.
