# Independent scientific semantic review of verifier ownership refactor

Reviewed 2026-11-28 against authoritative `git show HEAD:crates/core/src/verification.rs`, with the current `planning/plan.rs` and `verification/{engine,evidence,pcr,terminal,sources,projection,output}.rs`. Also compared the protocol PDF entry points/rendering body against `HEAD:crates/core/src/pdf.rs`.

**Result: no scientific semantic regression found in the inspected refactor.** No dropped checks, altered execution ordering, changed budgets, lost partial results, broadened physical authority, changed plan-revision preimage, or imported-report protocol construction was found. This is a source review; no compiler or tests were run for it.

## Comparison method and coverage

The original verifier was normalized with the current canonical namespace rewrite (`ns.transform(original, is_core=True)`) to separate module-path changes from code behavior. Comparisons used lexical tokens that retain string/character literal contents; whitespace outside literals and comments were ignored. Every residual token delta was inspected rather than treating all syntactic changes as equivalent automatically.

All 49 original top-level functions are accounted for:

- 32 retain exact complete function tokens after namespace normalization.
- Nine have exact bodies with only restricted visibility changes: evidence_bytes, retain_topo_source_incompleteness, validate_topo_shape, region_map, reverse_reference_map, project_into, project_into_with_collapse, carry_unmapped_evidence, and tailed_document.
- evaluate_captured now owns its state through Execution; preflight/report initialization and post-reconstruction status/error handling compare exactly after the intended ownership change.
- export_protocol_html now invokes the fresh context constructor and its HTML method; the corresponding constructor/renderer bodies were independently compared.
- Six former functions are now owner methods: operation_sources, operation_sources_mut, plan_revision, reconstruct, checked_protocol, and export_protocol_html_checked.

All 27 original struct/enum declarations, including their serde/schema attributes and field/variant order, compare exactly modulo restricted visibility. LigationOutcomeEvidence, MaterialStateEvidence, DerivedOriginEvidence, EvidenceCounter, and ProtocolAppendix validation implementations retain exact tokens after namespace normalization.

## Execution state, error prefixes and budgets

Execution owns VerificationReport and its physical-prefix vector while borrowing the existing final Golden Gate capture and CleavageBudget mutably. Its reconstruct method destructures mutable references to the same report products/findings/step/evidence vectors. It does not clone a budget, substitute an imported report, reconstruct from target DNA, or create a separate product representation.

The reconstruct body from target.validate onward retains all original substantive tokens and their order after changing operation_sources(step) to step.sources(). Residual differences are formatting-only: trailing commas and expression block wrappers around existing closure/match-arm expressions. No check or computational expression was removed or replaced.

Specifically retained:

- Target validation and separately bound target revision; step/input limits; valid unique input IDs; exact input revision checks.
- CheckedTopo final-step restriction and reference-only ancestry refusal for physical PCR/annealing/TOPO/homology/resection substrates.
- Forward/cyclic dependency rejection, original-input unknown duplex initialization, and refusal to invent missing tracked physical state from step metadata.
- The complete operation match, its bounded evidence collection and every per-operation guard.
- Reference-only flag propagation, refusal to attach physical state to a reference-only result, and restrictive portable reference markers.
- Tracked modification summary retention; 32 million base materialization accounting including both physical strands; 64 MiB product-plus-evidence accounting; the existing 32 MiB step-evidence bounds.
- Plan/engine metadata insertion and the final publication order: product push, physical push, then reference-only-prefix push, after the same materialization checks.

In evaluate_captured, preflight still happens before plan hashing. Invalid preflight returns empty products/physical state and unbound revision strings. Successful reconstruction still compares the final product to the independent target and uses ConsistentUnderAssumptions versus Incomplete according to the same incomplete findings. Unsupported failures still become Incomplete; other errors retain Invalid. Execution is destructured after either result so the actual completed product/physical prefix is retained, exactly as the old borrowed vectors were.

verify_with_physical, verify_with_tracked, verify_with_golden_gate_trace, verify_with_tracked_budget and evaluate retain exact function tokens after namespace normalization, including mutable capture/budget sharing.

## Plan identity and source visitors

Operation::sources and sources_mut have exact old visitor bodies after replacing the parameter identifier with self. Variant coverage, source order, oriented-source unwrapping, and mutable-reference behavior are unchanged.

DesignPlan::revision has exactly the previous plan_revision body after the identifier-only substitution plan -> self. Its preimage remains serde_json::to_vec(self), with the same blake3 hash and hex encoding. The seven plan/source/operation type declarations and their derives/serde annotations are unchanged, including field and enum variant order. No internal version or semantic serialization change was introduced by moving the types.

## Protocol authority and output

CheckedProtocol retains the exact same private fields and has no deserialization derive. Its new constructor has the exact old checked_protocol body: borrow the supplied plan and target, call fresh verify(plan, target), and initialize no appendix. There is no constructor parameter or setter for an externally supplied VerificationReport. with_appendix still validates the documentary appendix before attaching it; plan/target/report getters remain immutable.

The HTML renderer method retains all old substantive renderer tokens after changing checked to self. The only remaining token deltas wrap two existing TOPO-contour format expressions in match-arm blocks. Status text, warnings, incomplete results, processing/reference evidence and HTML escaping remain the same.

PDF preflight, export_protocol_pdf, and export_protocol_pdf_checked bodies compare exactly after canonical namespace/constructor rewriting. export_protocol_pdf_checked_inner differs only by two inserted trailing commas. It still performs preflight, uses checked.report(), initializes layout from that fresh report status, and preserves draft/evidence rendering. Other current constructor uses in quantities and combinatorial export use CheckedProtocol::new, and no alternate literal construction outside verification/output.rs was found.

## Review artifacts

- Original authoritative snapshot: `/tmp/collagene2-verification-head.rs`.
- Canonicalized original: `/tmp/collagene2-verification-head-canonical.rs`.
- Function comparison inventory: `/tmp/collagene2-verifier-review-comparison.json`.
- Read-only comparison scripts: `/tmp/collagene2-review-verifier.py`, `/tmp/collagene2-verifier-token-deltas.py`, `/tmp/collagene2-verifier-format-deltas.py`, and `/tmp/collagene2-verifier-pdf-review.py`.

The namespace manifest was used as comparison input only. No source files, repository metadata, runtime data, or protected documentation were modified or read beyond the assigned source/mapping scope. The integration owner's compiler and test results remain the execution evidence for the complete batch.
