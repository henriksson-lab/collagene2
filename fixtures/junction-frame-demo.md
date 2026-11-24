# Junction DNA and explicit reading frames

This short synthetic example checks sequence/codon arithmetic, not a feasible
wet-lab assembly. Open `junction-frame-demo.fa`, activate `junction-target`, then
choose **Assemble** and the overlap workflow. Use two whole-linear source rows:
`junction-left` then `junction-right`, both forward. Choose a linear product and
one overlap of **1 base**. This deliberately lies outside ordinary method length
guidance; that guidance must remain visible and advisory.

The shared homology is the single `T` at product interval `[3,4)`. The expected
product is `AAATGATTT`, with only one copy of that shared base. Review and verify
the candidate, then explicitly request junction inspection in the plan panel.

With a flank long enough to include the whole nine-base product, forward
**window-relative frame 0** reads `AAA TGA TTT`. Genetic code 1 produces `K*F`;
code 2 produces `KWF`. The junction boundary after the shared interval lies
inside `TGA`. Changing the requested genetic code must clear the old inspection
before a fresh result is shown. These are translations of actual window DNA,
not inferred gene intent, expression or fusion function.

Frames are relative to the oriented displayed DNA window. Their offsets are not
silently presented as product-origin frames. Original product ranges, the shared
homology interval and the boundary offset remain explicit. Separate CDS evidence
uses the actual retained annotation geometry and phase when supported; reviewed
or otherwise uncertain annotations do not become authoritative frame assertions.

The JSON download retains the exact plan/target/options context and independent
verification outcome. A wrong target may still have an inspectable reconstructed
product; failed final reconstruction must not substitute a preceding preparation
as that product.
