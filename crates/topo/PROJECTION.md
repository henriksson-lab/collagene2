# Source-bound paired-contour projection

`collagene_topo::projection` is a view proof over a freshly reconstructed physical
TOPO graph. It is not an additional reaction, repair operation, or kit prediction.
The original transfer engine remains unchanged.

Call `project(substrates, transfers, engine_limits, selection, options)`. Selection
pins the original engine batch revision, component number and an immutable source
`Call`. The projection revision separately binds the exact sources, transfers,
engine limits, seed, component, projection policies and projection limits. No target
sequence participates. A circular reference begins at the seed, in its actual
5′→3′ sense. A linear reference begins at its physical endpoint; the seed position
is reported. Changing source roles or chemistry invalidates selection.

For each paired call, forward continuation must agree between its actual covalent
successor and the paired counterpart of its partner's covalent predecessor. The
inverse condition must hold as well. The proof requires one antiparallel path or
cycle covering every basepair of the selected component exactly once. A missing
bond may be bridged geometrically only by the opposite strand's actual covalent
bond; the missing bond remains a nick. Branches, cross-track backbone connections,
inconsistent cycles, ambiguous DNA and unpaired calls withhold the view.

A proven `DuplexView` contains both actual strand sequences, compressed source-call
maps, source-fixed origin, external end chemistry and every nick's two source calls
and terminal chemistry. Bottom DNA comes from observed source partners, rather
than filling in complements. `Phosphotyrosyl`, blocked and unknown termini remain
explicit. Circular contour topology does not imply covalently closed DNA. This
view cannot be silently substituted for a simpler duplex type that drops nick
chemistry, flaps or activations.

Every result retains the complete fresh physical report, including all unselected
components and byproducts. Selecting a fully paired component does not delete a
separate free strand. Selecting a component containing a directional flap withholds
projection and leaves the flap graph intact. Ordinary TA nicked circular products
can be projected without sealing either nick.

Defaults permit geometric inspection of retained state when a requested transfer
was not applied, and retain unknown chemistry. The two explicit strict policies
can instead withhold such views. Neither default proves successful cloning.

The original engine's source/work/output limits apply first. Projection additionally
bounds pairs, conservative work cells and scratch bytes before source-sized proof
allocation, mapping segments, and the complete serialized result. A budget outcome
is withheld evidence, not absence of another geometry. Defaults are 50,000 pairs,
16 million work cells, 32 MiB scratch budget, 8,192 mapping segments and 32 MiB
whole-report output; hard limits are checked before hashing.

Serialized reports are unverified claims. `replay` reruns both source-based transfer
and contour proof and compares the complete result. A report alone does not restore
physical authority. The reference view does not model supercoiling, stability of
short duplexes, binding efficiency, host repair, strand methylation or proprietary
kit preparation. Primary chemistry references remain in `data/primary_sources.json`.
