# Ligation alternatives in protocol output

The normal protocol PDF renderer handles `Ligate` and the fresh
verifier's `LigationOutcomes` evidence. It does not run a second target-selected
alternative search: the protocol's ordinary fresh native reconstruction supplies
the actual intended material and the bounded documentary alternative report.

Every retained product, including the separately named intended path, is printed
with its complete ordered/oriented source path, reference length/topology, actual
top and bottom 5′→3′ strand sequences, bottom-reference offset, free-end chemistry
where applicable, retained nicks and modification-state identity/counts/generation
provenance. Source operand keys retain exact source and actual-material revisions.
The intended target remains independent of the enumeration. Alternatives do not
prove target success, isolation or reagent activity, and no probabilities appear.

The text states the exact declared scope: one copy of every operand, the requested
topology and nick policy, with circular coordinate anchoring and equivalent linear
views described in the native report notes. Partial assemblies, other copy counts,
vector-only products, multimers, concentrations, kinetics and yields are outside
this scope. It prints requested budgets, actual checks/work/retained bases,
incompatible/unresolved/withheld counts, completion within this scope, every
reached limit and all retained warnings/unresolved examples. The native bounded
unresolved-example excerpt remains explicitly labelled with its total/truncation;
it is never described as an exhaustive list of unresolved possibilities.

Each ligation step also embeds two ordinary PDF file attachments:

- `ligation-step-NNNN.json`: the complete typed evidence, including bindings and
  every retained intended/alternative product and report limit.
- `ligation-step-NNNN.fasta`: every retained top and bottom strand, wrapped at80
  symbols and labelled by step, intended/alternative index, topology and strand
  orientation. FASTA is sequence-only: end chemistry, nicks, modification state and
  the full declared scope remain in JSON and printed text. The bottom strand is
  its own5′→3′ sequence, not a second forward-reference sequence.

Attachments use pdf-writer's embedded-file/name-tree APIs; no filesystem,
subprocess, new dependency or active PDF script is involved. Viewers that support
PDF attachments can extract them. The ordinary reconstruction JSON export remains
available independently of viewer attachment support.

Existing limits remain256 pages,512KiB printable text and16MiB total PDF, including
attachments/font/content/object overhead. Writers bound evidence bytes before
allocation grows beyond the file allowance; the final assembled PDF is checked
again. No partial PDF or truncated alternative list is returned. An oversized
ligation export explicitly directs the user to full reconstruction JSON. This is
an explicit bounded refusal, not permission to silently omit a long sequence or a
late alternative. Embedded evidence is documentary and never restores native
physical-material authority when reimported elsewhere. This is not PDF/A.

The compact web processing panel shows exact scope/status and lets the user inspect
any retained product, including the intended path. The selector changes only the
inspection. Complete strand sequences, geometry/state, all bounds and source pins
are available in nested details; no product-acceptance action or target-success
claim is added to this panel.

`core/tests/ligation_outcomes_pdf.rs` authors independent PDF text/CMap/xref and
length-delimited attachment readers. Literal PCR→phosphorylation→blunt circular
ligation checks both the independently specified intended sequence
`AGTCCGGACTGA` and reversed-insert alternative `AGTCCGTCAGTC`, and byte-equivalent
complete JSON/FASTA evidence. A separate bounded-product case retains actual nicks
and states incomplete enumeration/withheld products. Tests are run by the parent;
this source worker did not execute a compiler, browser or application request.
