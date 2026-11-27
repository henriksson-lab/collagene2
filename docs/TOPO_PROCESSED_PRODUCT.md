# Source-seeded TOPO contour and faithful material bridge

All transfers, flap cleavage and nick sealing use the one canonical ordered
transcript and fresh private cursor described in [physical actions](TOPO_PROCESSING_ACTIONS.md).
`physical::projection` proves a selected paired contour over the final actual
graph. Its `Selection { source_batch_revision, component, reference_seed }` pins
the canonical physical revision and original source call. The seed chooses strand
sense and circular origin; intended target DNA never guides the path.

`source_batch_revision` is `TopoReport.physical.batch_revision`, distinct from the
outer source-adapter `TopoReport.batch_revision`. The whole physical report,
including other components, free calls, conditions, outcomes and activation
history, stays in `ProjectionReport.physical`. Strict requested-step completion
and strict known-terminal chemistry are separate options. An unpaired directional
flap remains withheld until an explicitly supported cleavage actually removes its
bond from the selected contour; its excised calls remain elsewhere in the graph.

There are two different outputs:

- `core::topo::Product` supplies a restrictive reference molecule with optional
  annotation correspondence. This does not establish reusable physical material.
- `core::topo_actions::Report` is the sole faithful contour-to-duplex bridge. Its
  `Request` contains `preparation: topo::Request`, canonical physical `selection`
  and `options`. It has no separate processing or sealing phase request.

The bridge report retains `engine_version`, `context_revision`,
`preparation_revision`, `request`, `input_sources`, `sources`, `preparation`, canonical `projection`,
`duplex`, `tracked`, `product`, `annotations` and `bridge_notes`. There is no second
processed-product API or processed contour representation. Root `input_sources` bind the external caller identities; expanded `sources` bind the native graph/annotation namespace. A graph-only continuation legitimately has empty roots but retained expanded dependencies; both identities are checked independently.

## Representability and authority

Only a complete aligned paired contour can become an ordinary duplex. Linear
termini retain exact representable states. Internal ordinary nicks require
explicit 3′OH/5′P chemistry. Retained phosphotyrosyl, 3′phosphate, blocked 5′ or
unrepresentable nick chemistry stays in the graph with an explicit bridge refusal;
it is never normalized to a convenient end. Circular placeholder termini are
unused; actual backbone closure comes from the graph.

The bridge retains Unknown internal modification calls. Annotation correspondence
asserts neither chemical inheritance nor function. The standalone product molecule
carries the physical-evidence restriction and remains a reference when exported or
reimported; serialized `tracked` data and positive report JSON do not establish
native material. A verifier must freshly execute the source-bound transcript.

Component selection is an explicit view, not experimental isolation. Other
components and released/excised material remain in the full report and in native
[graph records](TOPO_GRAPH_HISTORY.md). Supported saved graphs can be continued
without flattening them through this narrower duplex bridge.
