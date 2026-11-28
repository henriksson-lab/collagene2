# Construct QC: initiation context, ORFs and causal junctions

Source review and implementation batch: 2026-11-28. No compiler, test, browser or
application HTTP checks were run by this worker. The new independent fixtures are
source contributions to the parent's combined native/WASM verification batch.

## Implemented in the existing owner

`crates/core/src/construct_qc` remains the single source-bound QC owner. Its current
engine is `construct-qc-3`, options schema is 3, and older internal options are
refused. There is no additional plan, verifier or crate. Options, captured source,
annotations, genetic-code tables, profiles and declarations enter the existing
request revision and fresh `Report::validate_for` replay. Neither new check adds an
overall biological pass.

| Check | Input and evidence | Interpretation |
| --- | --- | --- |
| Initiation pattern | `initiation_context::ContextClaim`, exact `PatternProfile`, `ContextEvidence` with observed oriented sequence, mapped coordinates, primary citation and verdict | Literal conformity to a declared host/transcript-context pattern. A mismatch fails the **pattern predicate**, not expression. |
| Unannotated ORF candidates | `unannotated_orfs::ScanOptions`, canonical ORF enumeration, `ScanEvidence`/`Candidate` with location, boundaries and annotation relation | Informational candidates, never errors or evidence of unintended expression. Unresolved annotation coverage is indeterminate. |

`Options.initiation_contexts` is a required bounded vector; `Options.orf_scan` is a
required nullable field. Reports retain corresponding evidence. No missing-field
fallback or superseded schema decoder is introduced.

### Initiation context

`PatternProfile` captures id, content revision, declared host context, genetic
code, IUPAC sequence pattern, the zero-based A-of-AUG index and documentary
citation. `compute_revision` pins every scientific/applicability field;
`validate` checks that pin. Citation URLs identify provenance; an arbitrary URL or
caller declaration does not authenticate a biological reference.

The provided `kozak_1986_pattern()` captures the literal DNA sequence `ACCATGG`
with AUG starting at offset 3, for a declared mammalian context and code 1. This
sequence comes from the primary preproinsulin mutagenesis study, not from an
expression classifier. The published study is context-specific; conformity is
neither a necessary nor sufficient expression condition. [Kozak, Cell 44:283–292
(1986)](https://doi.org/10.1016/0092-8674(86)90762-2). The later 699-vertebrate-mRNA
survey supplies broader documentary context, not a universal host model.
[Kozak, NAR 15:8125–8148 (1987)](https://doi.org/10.1093/nar/15.20.8125).

The rule requires a current declared, complete CDS and an exact single interval,
optionally complemented. Existing CDS translation establishes the boundary and
code. The context window may wrap a circular source once, uses biological
orientation and accepts T/U mask equivalence. Every source expansion contained in
the pattern means conformity; disjoint masks mean mismatch; overlapping masks
with unresolved expansions mean indeterminate. A non-definite AUG, including a
recognized non-AUG initiator in another code, is inapplicable/unresolved for this
profile. Missing linear flanks, a window longer than one source turn, partial,
unreviewed, fuzzy or compound CDS paths, and a code mismatch remain indeterminate.

Current source DNA around a CDS does not prove mature transcript boundaries,
splicing, promoter use, transcript existence or host applicability. These are
explicit declared assumptions, displayed in the evidence. Compound mature
transcript context needs an explicit transcript-coordinate owner before support;
it must not borrow adjacent genomic DNA across an intron or infer a UTR.

Admission is at most 64 claims, 128 bases per pattern/window, bounded ids,
applicability/rationale and citation fields, existing 2 MiB options and report
caps, and a 2M-base current source. Evidence and citation copies charge the shared
QC report budget. No whole-genome context copy or sequence alignment is used.

### ORF candidates

The rule calls `sequence::translation::find_orfs_validated`, retaining its current
start policy, partial policy, code, strand, ambiguous-stop and circular one-turn
semantics. It adds no codon decoder or matcher. Users declare a positive amino-acid
threshold and 1–128 candidate cap. The underlying owner scans linearly one strand
at a time and has an 8M-amino-acid aggregate materialization bound. Proteins are
not duplicated in QC evidence. Annotation comparisons are bounded by 128
candidates and 256 source features; output remains within the shared report cap.

`enumeration_complete` means the canonical enumeration completed **under its
specified policy**. It is false on a result/materialization cap. Ambiguous
possible starts are not enumerated by the definite-start owner; even a completed
scan is not exhaustive detection of all potential coding content. Ambiguity and
context-dependent stop tables are separately visible and give an indeterminate
summary. Available candidates remain visible after truncation. Partial and
uncertain-stop candidates retain their explicit boundary state.

Coverage requires a current declared complete CDS in the same selected code,
containing the entire candidate on the same strand and codon register. Identical
exact biological paths also cover origin-spanning candidates. Nested accepted
starts are retained. Simple interval overlap, feature names and genomic
`start % 3` are insufficient. Potentially relevant unreviewed or unsupported
paths make coverage unresolved; they cannot establish an unannotated result.
`NoDeclaredCdsCoverage` is a literal relation to current annotations, not proof of
biological unintendedness, gene absence or expression. No candidate gets a
biological fail or a functional pass.

Independent literal fixtures cover forward/reverse/circular windows, RNA,
IUPAC subset/disjoint uncertainty, missing flank, partial/unreviewed CDS,
alternative start and code mismatch; nested ORFs, wrong-register overlap,
review-required coverage, reverse/origin ORFs, result caps, ambiguous starts and
stops, one-turn partials, source/options/annotation changes, tampered content and
required-null/current-version admission. Fixtures use known literal expected
coordinates and sequences; their implementation is not the expected oracle.

## Remaining causal-junction owner work

Existing `JunctionClaim` establishes only a current declared boundary/register.
Its history ids are documentary. `documents::workspace::RevisionEvent` retains
source revision, ordered parents and a descriptive operation string, explicitly
not executable replay authority. Event ids alone cannot prove that an operation
introduced a frameshift or stop.

The available stronger foundation is `material::lineage::RunRecord`: the current
`DesignPlan`, exact target/input/product `EventBinding`s, step indices and bounded
workspace attachment. `load_run` and fresh `replay_run` can recover and verify this
record against immutable source events. Recorded assessments are historical
statements and must not substitute for fresh verification.
`cloning::assembly::junctions::inspect_plan` already reconstructs the expected
product and its junction geometry through the current verification engine,
including method-specific prepared fragments and CDS windows. Reuse this owner;
do not implement a second assembly simulator in QC.

A coherent next adapter should:

1. Resolve one selected saved run/product event from a validated workspace;
   retain exact attachment, run, event, step and source bindings in the report.
2. Freshly replay the plan and inspect its junctions, require the chosen product
   to equal the current molecule under the actual source identity owner, and
   explicitly reject another branch or a changed source/annotation/options pin.
   Do not compare unlike revision encodings: the junction owner's bounded-JSON
   product hash and `Molecule::revision()` are different identity domains.
3. Establish a specific parent/preparation-to-product **coding-path mapping**
   before claiming causality. Reused labels, overlapping features, equal sequence
   or inferred annotation roles do not supply this mapping. Use typed operation
   coordinates/preparation provenance where retained; ordinary descriptive edit
   history is insufficient.
4. Compare the corresponding upstream and downstream codon registers and codons
   under the same pinned translation model. A mapped one-base insertion can
   establish a register change; a mapped three-base insertion preserves register
   but may still alter a codon/stop. A pre-existing stop or already-disrupted
   parent must not be attributed to the new junction.
5. Accumulate current structural findings independently when ancestry, feature
   correspondence, an unsupported operation or ambiguous codon blocks causal
   interpretation. Label the causal rule indeterminate, not absent or passing.

This adapter should be a cohesive subowner in current QC/material review, not a
new general plan representation. Keep one selected bounded replay per request,
reuse the current junction owner limits (including 128 junctions), admit aggregate
source/report bytes before clones, and avoid replaying arbitrary entire history
DAGs. General edit causality needs retained typed coordinate transforms where
none exist today; operation text must never be parsed as trusted instructions.

Required independent fixtures for that adapter: mapped +1 versus +3 insertions;
a reverse CDS; an origin-spanning coding path; pre-existing versus newly crossing
stop; two simultaneous changes without a unique attribution; changed source
revisions; selected sibling products with equal DNA but different event identity;
truncated junction inspection; missing/unsupported operation mapping; altered
attachment hash and restored stale report. None should assert experimental
assembly, expression or a causal claim from one current boundary alone.


Parent integration checkpoint: the shared native/WASM batch and 349 unique focused
cases pass. See [source and execution evidence](cloning-extensions-current-integration.json)
for the exact scope, initial catalogue-pin corrections and final metadata-only
wording edit. No graphical or application HTTP acceptance is claimed.
