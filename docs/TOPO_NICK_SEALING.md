# Explicit downstream nick sealing

TOPO strand transfer and subsequent ligase treatment are separate source-backed
steps. `collagene_topo::nick_sealing` freshly reconstructs the original transfer
reaction, then evaluates only explicitly selected pairs of source calls. It does
not accept a serialized product graph as starting material.

The declared model requires a free 3′ hydroxyl and 5′ phosphate at actual exposed
chain ends, exact paired DNA and a continuous antiparallel template across the
nick. These terminal requirements follow the documented
[T4 DNA ligase reaction](https://www.neb.com/en/products/m0202-t4-dna-ligase).
The model deliberately covers an ideal templated nick under an explicit ligase/
conditions assumption. It does not claim all ligases behave identically or rule
out reactions outside the supported geometry.

Successful treatment adds the selected backbone edge. Every original nucleotide,
source placement, pairing, released activation and displaced call remains
accounted for. A selected step with unsupported geometry, blocked chemistry or
unknown state preserves its incoming graph and records its outcome. Later steps
operate on that explicitly retained branch; it is not an exhaustive reaction
mixture. Duplicate or already consumed ends cannot be silently joined twice.
No kinase, phosphatase, new base synthesis, flap removal or host repair is inferred.

The original physical transfer report is retained alongside the treated product.
Two explicitly sealed TA nicks can close both backbones under the declared model;
a directional displaced flap is still present after treatment. Neither result
implicitly becomes an ordinary duplex substrate or establishes agreement with an
independently intended target.

## Application and MCP

`core::topo_nick_sealing::Request` contains `preparation` (the original TOPO recipe)
and `sealing` (the explicit ligase assumption and selected nicks). Source molecule
identities, preparation, transfer and treatment choices bind the context revision.
Preparation, including exact PCR if requested, runs once. Fresh replay reconstructs
the source graph and treatment; imported evidence never restores authority.

The existing TOPO orientation review opens a downstream nick-treatment panel.
The user loads an inspected orientation's actual nicks, selects which ones to
seal, enters the assumption and requests review. Changes to source events,
orientation or treatment invalidate pending work and displayed evidence. Results
retain all components, remaining nicks, unpaired runs and terminal states; full
JSON evidence can be downloaded explicitly.

Worker63 appends job147 `SealTopoNicks` and result139 `TopoNicksSealed` with ordinary
message limits and historical-version rejection. MCP161 adds `seal_topo_nicks`
using immutable source bindings and the same typed request. This operation returns
a graph report, not a newly verified document or an accepted cloning candidate.
Existing plan/candidate/science versions and TOPO layouts remain unchanged.

Remaining F15 work includes accounted directional flap processing, fresh contour
proof and tracked-duplex composition after processing, durable graph-only workflow
handoff, broader method profiles and final delivery/GUI acceptance. Adding a
ligase step does not imply these requirements are complete.
