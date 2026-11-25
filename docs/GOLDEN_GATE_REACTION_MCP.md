# Whole-reaction Golden Gate MCP workflow

`inspect_golden_gate_reaction {inputs,target,request}` takes exact original
bindings and an independently declared target. Each row explicitly chooses
`WholeSource` or `Prepared {preparation,treatment}`. Complete ExactPcr material
includes the entire amplicon; a selected catalog digest remains intentionally
isolated material. The response reports every compact candidate summary,
cut/carryover scope, caps, source-role membership and independent target matches.
The full artifact `golden-gate-reaction.json` contains the versioned evidence kind,
exact bindings, request, material map, target, probe plan and actual native pool.

`probe_plan_revision` identifies diagnostics for native index0 execution, which
is not automatically registered as a candidate. `probe_verification` must not be
read as a chosen reaction result. Recutting and unavailable probe index0 still
retain the actual pool; native preparation failure does not invent one.

`select_golden_gate_reaction_candidate {artifact_id,selection}` requires the exact
context/index/directed+physical keys from that artifact. It checks the retained
source and target references, then atomically saves an unverified plan and binary
`candidate.cg2plan`. `artifact_id` and `candidate_artifact_id` identify the same
bytes. Use `verify_plan` on the new plan revision before product review, capture
or HTML/PDF instructions. Native reconstruction, not serialized pool state,
supplies physical material. Complete cleavage and finite one-copy assumptions
remain explicit, with no experimental or exhaustive-mixture verdict.

Independent fixtures declare a40-base target before reaction analysis. A circular
vector contributes both its backbone and stuffer-side digest piece; a complete
PCR amplicon contributes three pieces; a third uncut circle is separate carryover.
The six-fragment pool retains self-products, recuttable backgrounds and the full
intended construct. Selection and fresh verification distinguish that construct
from self-products. One-state/minimum-two-fragment caps retain an empty but actual
incomplete pool; source failures and changed identity keys cannot nominate DNA.
Validation results are in [implementation status](IMPLEMENTATION_STATUS.md).
