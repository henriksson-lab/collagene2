# Ordered TOPO transfers and physical processing

`physical::Step` is the current chronological representation:

- `Transfer(Transfer { donor: Call, acceptor: Call, geometry })` consumes an actual
  exposed phosphotyrosyl donor and an actual exposed 5′OH acceptor.
- `Process(processing::Action::CleaveThreePrimeFlap { ... })` declares one exact
  supported hydrolytic cut.
- `Process(processing::Action::SealNick { ... })` declares one exact templated
  3′OH/5′P ligation.

`Call { strand, index }` addresses an immutable source nucleotide. The executor
uses its current private backbone, pairing and effective terminal chemistry.
Donor motif checks follow the actual upstream backbone, rather than assuming a
source-string end or subtracting reference coordinates. `StepEvidence` mirrors
the ordered declaration, using global transcript indices. A flap's
`displacement_step` identifies the earlier actual transfer in the expanded source
and current transcript, including processing positions.

## Declared chemistry

`IdealDuplexAnchoredThreePrimeHydrolysisV1` is an explicit ideal hydrolysis
hypothesis. It is not a named nuclease, kit protocol, cellular-repair model or
experimental-success prediction. The selected bond must belong to the exact
recorded duplex-anchored terminal four-call flap. Applied cleavage creates an
upstream 3′OH and downstream 5′phosphate. The detached calls, their distal original
chemistry and their own component remain present. A target never selects the cut.
Unsupported or blocked actions preserve the incoming state and explain why.

The new 3′OH is not a TOPO activation; the new 5′phosphate is not 5′OH. Later
transfers and nick treatments check these actual states. Consumed activation
identities remain consumed even when other backbone bonds are subsequently cut.
No phosphorylation, dephosphorylation, fill-in, reactivation or host rescue occurs
without a separately supported explicit operation.

## One fresh execution

The private cursor holds original-call backbone/pairing/displacement provenance,
terminal overrides and consumed/released activation state. Transfers and
processing mutate that same cursor in order. A shared work ledger also covers
orientation branches, source checkpoints and fresh contour projection. Source,
step, work, memory, output and report reservations precede their corresponding
allocations or operations. There is no report-to-cursor reconstruction API.

`Report.request` retains the canonical declared steps and conditions;
`Report.steps` retains their typed evidence. `all_requested_steps_applied` covers
both transfer and processing. Fresh output independently audits every original
call exactly once, including split source intervals and free flaps. Readonly
binding validators check bounded correspondence; they cannot manufacture fresh
material or replace scientific replay.

## Source continuation and products

`core::topo::Request.steps` and native graph records use the same transcript.
A saved graph executes its whole retained prefix before the newly declared suffix.
Checkpoints use `after_steps` and `step_map` and hash effective end chemistry,
sealed backbones, all components and displaced-call ownership.

The sole `physical::projection` proves a source-seeded contour over this actual
fresh graph. `core::topo::Product` is a reference view. The distinct
`core::topo_actions` bridge can produce private tracked duplex material only when
all relevant strand/end/nick chemistry fits that model faithfully. The entire
unprojected graph remains in the evidence; selection is not physical isolation or
a reaction-mixture winner. See [native records](TOPO_GRAPH_HISTORY.md),
[faithful material bridge](TOPO_PROCESSED_PRODUCT.md) and
[plan reuse](TOPO_PROCESSED_PLAN.md).

The independent fixture in `crates/core/tests/cloning/topo/graph/topo_graph_continuation.rs` specifies a
52-call directional graph: first transfer, actual flap hydrolysis and one nick
seal are captured; after restart the remaining original activation is transferred
and the other nick sealed. The result retains 48 paired circular calls and the
free four-call component. Newly generated OH/P endpoints, duplicate activation
use, changed selections and duplicate native source roots are checked separately.
