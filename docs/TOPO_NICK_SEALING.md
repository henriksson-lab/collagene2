# Explicit nick sealing in the TOPO transcript

Nick sealing is a `physical::Step::Process(processing::Action::SealNick)` in the
same ordered transcript as phosphotyrosyl transfers and explicit flap cleavage.
There is one fresh private cursor and one canonical physical report. The
`nick_sealing` module contains the selected-call/evidence types and eligibility
rule; it is not a separate reaction executor or report format.

The action declares `NickSelection { three: Call, five: Call }` and a bounded,
explicit `ligase_assumption`. Both calls identify actual immutable source
nucleotides. Eligibility requires exposed 3′OH and 5′phosphate, exact paired DNA,
and a covalently continuous antiparallel template across those termini. These
terminal requirements follow the documented [T4 DNA ligase reaction](https://www.neb.com/en/products/m0202-t4-dna-ligase).
The supported profile is an ideal templated nick under declared conditions, not
an efficiency prediction or a claim about every ligase or unsupported substrate.

Applied sealing adds only the selected backbone edge. It does not alter pairing,
sequence, activation identities, other nicks or other components. Unknown
chemistry remains unknown. An activated 3′ donor is a distinct TOPO transfer;
sealing cannot consume its activation. A 5′OH does not become phosphate, and a
blocked end does not become hydroxyl. Repeated or already joined endpoints retain
the incoming graph and report the blocked outcome.

Hydrolysis-created internal ends use their actual per-call 3′OH/5′P overrides.
They do not inherit the original distal source terminus. Later actions see the
actual new backbone and chemistry, including already sealed bonds. No report,
projected molecule, metadata marker or intended target restores this state.

`StepEvidence::Process(ActionEvidence::NickSealing(...))` records the global
transcript index, selected calls, actual terminal states, opposite template calls,
outcome and reason. Operator conditions remain in `Report.request`; canonical
outcome summaries cover every requested step. Projection can inspect a retained
branch permissively, but an unapplied action remains unapplied and inherited
incompleteness is not erased by sequence agreement.

Use the canonical `topo::simulate` graph inspection and, when a selected contour
is faithfully representable, the sole `topo_actions` tracked-material bridge.
[Graph capture/replay](TOPO_GRAPH_HISTORY.md) preserves this same transcript across
restart. There is no separate sealing context, phase decoder or processed-plan
alias.
