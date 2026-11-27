# Source-backed TOPO adapter

`core::topo` reconstructs native source recipes from immutable
`Input { role, molecule }` records. Activation, terminal chemistry, pairing and
reaction conditions are explicit declarations, never inferred from a name,
vendor, annotation or imported physical-summary metadata.

A `Prepared` material describes stored 5′→3′ strands assembled from bounded
original-source intervals, with explicit pairing spans and termini. Intervals may
reverse-complement a reference region or cross one circular origin. References
are generation templates under a stated hypothesis, not proof of atom inheritance
or manufacture. Unrepresented source intervals remain documentary evidence; they
are not invented reaction waste.

`ExactPcr` declares actual ordered primers and optional subsequent single-base 3′
tailing. Native geometry must establish the product from the actual permitted
template. Supported saved Run/Derived templates are freshly reconstructed; their
modifications do not transfer to new synthesis. The current TOPO template adapter
requires complete aligned complementary strands without nicks. A reference-only
projection or unsupported staggered source cannot bypass this boundary.

`SavedGraph` is a distinct closed native graph source. It retains the whole graph,
source snapshots and chronological prefix rather than pretending that an activated
vector or intermediate is an ordinary duplex. Native ownership, exact events and
one aggregate ancestry ledger precede copying; fresh execution establishes its
state. See [graph history](TOPO_GRAPH_HISTORY.md).

## One chronological physical transcript

`Request.steps` is a sequence of `physical::Step::Transfer` and
`physical::Step::Process` declarations. Transfer endpoints are immutable
`Call { strand, index }` handles, not changing output-chain indexes. The same
private cursor applies actual phosphotyrosyl transfers, supported flap hydrolysis
and selected 3′OH/5′P nick sealing in order. There is no report-to-cursor restoration
or separate processing/sealing phase. All components, cut fragments, unpaired
calls and consumed/released activation identities remain accounted for.

`simulate` returns complete preparation maps and a canonical physical report.
Its `request` and typed `steps` evidence retain the same chronological transcript.
`all_requested_steps_applied` covers transfers and processing. Work/memory/report
limits remain explicit, and unsupported/unknown steps retain honest outcomes.

`inspect_orientations` requires an empty new suffix and explicit vector/insert
strand pairs. It evaluates the two declared orientations after the actual retained
prefix, returning both full graphs and candidate-specific pending steps. It
chooses no winner and does not claim exhaustive real-mixture chemistry. Pairing
indexes in a Prepared material are local; expansion records their actual native
layout, role aliases and source-step mappings.

## Reference view and independent target

`product` freshly proves a source-seeded paired contour. Its `topo::Selection`
pins the outer adapter batch revision, component and original-call seed. The
physical projection inside that product separately binds the physical transcript
revision. Target DNA never supplies path choice, orientation or reference origin.
An unpaired flap withholds the contour until an explicit supported action changes
its actual graph; the detached calls remain elsewhere.

The one canonical `Product` contains optional conservative annotation evidence
and a restrictive reference molecule. PCR-derived or disrupted features require
review; dropped/unmapped originals remain source evidence. Coordinate
correspondence does not grant inherited chemistry or biological function.
Annotations are retained with the relevant `StepReport.topo`, not a second
annotation-only product API or generic always-empty report collection.

`check_target` compares a separately supplied target with this freshly proven
contour. Rotation and reverse complement are explicit options. Nicks, unknown
terminal chemistry, unapplied steps and other components each have their own
acceptance allowance. `Operation::CheckedTopo` is final-only; accepted sequence
agreement under that policy is not isolation, reaction success or permission to
reuse the reference as physical stock. `Operation::Topo` also remains a reference
view.

## Faithful ordinary-material bridge

`Operation::TopoActions` is the sole contour-to-tracked-duplex bridge. Its request
contains `preparation`, physical `Selection { source_batch_revision, component,
reference_seed }` and projection `options`. The physical selection revision is
`TopoReport.physical.batch_revision`, not the outer adapter revision.

Only faithfully representable strand, end and nick chemistry enters ordinary
native duplex operations. Unsupported terminal species, ambiguous contours and
unmodeled ancestry stay diagnostic. Internal modification calls remain Unknown;
source annotations do not infer a modification-generation model. Full graph and
byproducts remain in `physical_processing_reports` as typed `TopoActions`
evidence. Exported/reimported molecule metadata cannot restore that private state.

The shared bounded writer, source/ancestry checks, fresh replay and readonly
correspondence validators have distinct roles. Validators never rerun a supplied
report into authority. Native history and HTML/PDF retain original hypotheses,
chronological outcomes, source maps, actual nicks, detached calls and limitations.
See [physical actions](TOPO_PROCESSING_ACTIONS.md), [nick eligibility](TOPO_NICK_SEALING.md)
and [plan reuse](TOPO_PROCESSED_PLAN.md).
