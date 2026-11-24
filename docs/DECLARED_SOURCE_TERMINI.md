# Explicit original-source terminal hypotheses

`Operation::DeclareSourceTermini {template, left, right}` adds an explicit operator
hypothesis for an **original intact linear duplex reference**. It is appended as
operation tag 16; `Preparation::DeclaredDuplex {left,right}` is appended as shared
preparation tag 3. Existing variant indices and serialized structure layouts are
unchanged. The scientific engine version is science13; worker/candidate/plan
version gates distinguish the new interpretation.

Each existing `chemistry::Terminus` specifies 5′ phosphate, hydroxyl or unknown,
and 3′ hydroxyl, blocked or unknown. Left/right use the existing duplex terminal
convention. These are assumptions, not measurements or assertions inferred from
file metadata. Unknown groups remain unknown. `WholeLinear` retains its existing
unknown free-end chemistry and never automatically upgrades to a declaration.

Only `Source::Input` is accepted. The actual input must be linear, unambiguous DNA
that satisfies the existing intact-reference constructor. Single/mixed-strand
sources, circles, detached physical intermediates, reference-only/TOPO/variant
claims and active imported physical recipe/summary markers are refused. Marker
presence is restrictive even if a producer changes a value to `false`; it cannot
restore authority. Exact generated physical/reference-only COMMENT markers in
retained source text are also restrictive. A declaration cannot turn a plan
`Source::Step` into a new intact original, seal internal nicks, discard flaps,
synthesize DNA or replace a missing physical graph.

Fresh reconstruction first constructs the guarded original-reference duplex with
unknown state, then applies only the explicitly declared terminal groups. A
`TrackedDuplex` is created with **unknown positional modifications**, without a
new unmodified-DNA or methylation-generation claim. Subsequent checked digestion
therefore retains the existing unknown-modification restrictions. This operation
does not compose with `DeclareSourceModifications` on a descendant; both original
source assumption interfaces retain their explicit scope rather than resetting
tracked state.

The declaration preserves the source DNA, annotations, name, raw GenBank and
citations byte-for-byte. It adds explicit assumption/source-revision evidence,
fresh duplex-summary evidence and the normal verifier bindings/modification
summary. Those published summaries are restrictive detached evidence, not a way
to reimport phosphorylation or construct private tracked material. Reconstructing
the typed operation from its original source is the supported path.

The shared compiler emits the declaration as a real preparation step; row mapping
points to that step and preserves source identity/orientation. Existing optional
ligation `TreatEnds` and digest repair remain separate operations. Declared
5′-phosphate/3′-OH blunt sources can form sealed junctions in the existing model;
a singly unphosphorylated side can leave a nick only under explicit permission.
Unknown terminal groups withhold chemical conclusions; blocked groups do not
become accessible OH ends. Neither-strand-sealable junctions remain invalid even
when nick retention is allowed. An exact target match is conditional on the
operator hypothesis and never establishes experimental success.

Borrowed source/metadata/location bounds run before cloning or source hashing.
The helper `verification::declared_source_termini_preflight` is shared by native
plan execution and preparation validation. Existing verifier materialization,
annotation, request/output and scientific-engine limits remain in force. HTML
protocols print both declared terminal groups and the operator-hypothesis scope.
Focused literal fixtures cover sealed/nicked/unknown/blocked chemistry, retained
citations and unknown methylation, descendants/circles/forged marker refusal and
append-only operation/preparation tags.

A [small declared-ends ligation example](../fixtures/declared-ends-demo.md) independently specifies the source terminal choices and eight-base expected product.
