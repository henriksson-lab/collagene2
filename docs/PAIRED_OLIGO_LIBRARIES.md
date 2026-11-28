# Paired-oligo source libraries

A library choice declares its primary source event and an explicitly present,
nullable `anneal_bottom` event. An annealed choice requires both sources and an
explicit annealing candidate. Both source events, sequences and annealing options
participate in selection invalidation. Changing an equal-DNA source to a different
document event still changes provenance.

Selected constructs use one canonical input order: primary sources in bin order,
then bottom oligos in paired-bin order. Twelve target contributions can therefore
require 24 physical inputs. Equal sequences, or the same event used in two roles,
retain distinct input identities. Revision-graph parent nodes remain unique;
ordered physical role bindings retain multiplicity.

The bottom oligo is a preparation operand, not another intended target fragment.
Targets remain independently declared original regions or retained target events.
Fresh native annealing and assembly determine the product and compare it with
that target. Neither an imported verdict nor target DNA selects a reaction
candidate. Existing end-chemistry and source-admission rules apply to both oligos.

Choice preparations are reusable templates. Their local bottom-input label is
rebound to the selected construct's canonical role identity; the exact bottom
event is retained separately. Moving a choice between rows does not discard an
explicit candidate when its ordered source pair and options remain unchanged.
Changing either source, options or strand roles invalidates the candidate.
Incomplete catalogue choices may be saved, but cannot silently ignore a declared
bottom source by compiling another preparation.

Capture, archive restart, fresh replay, release, native run records and HTML/PDF
protocol output retain both ordered operands. Explicit plate publication includes
every source in the retained batch design, including unselected bottom choices.
Ordinary selected-plasmid publication continues to omit unrelated plate/batch
attachments. Missing dependencies cause rejection rather than DNA-based recovery.

Current contracts are source design 4, immutable batch kind
`combinatorial-batch-record-v3`/schema 3/context contract 4, release kind
`combinatorial-release-v2`/schema 2, active draft 20, library history 3, worker
schema 123 and sharing closure engine 5. Superseded internal encodings are
rejected; no migration or fallback decoder is provided.

The scientific, web-state, sharing and public MCP fixtures cover ordered operands,
equal-DNA distinct events, repeated-event roles, source-change invalidation,
explicit candidate selection, archive restart, release and protocol output.
Compilation/execution results belong in the [delivery audit](DELIVERY_AUDIT.md);
fixture source alone does not prove the current batch passes.
