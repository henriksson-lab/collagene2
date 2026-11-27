# Sequence similarity across selected documents

The native `sequence_similarity` module supplies the scientific core for F33
collection similarity search. Worker, web and MCP integration remain required;
this is not yet an application control.

Each query/subject has a bounded stable source key and its full typed molecule.
The report context binds those keys, exact molecule revisions, source order,
options and engine version. Equal sequence content in different documents is not
silently merged. Querying the same source included among the subjects is valid.
Subject keys must be unique within the requested collection.

## Declared search profile

The query is its **stored linear sequence**, optionally also reverse complemented.
Subjects are DNA or RNA of the same alphabet. Circular subject windows are rotated
under an explicit origin cap; each window contains at most one reference lap.
A circular query's partial matches crossing its stored origin are not separately
searched. Whole-molecule circular equivalence remains the separate comparison
workflow. Protein queries are unsupported in this profile; incompatible subjects
remain visible as unsearched rows.

The search returns one best observed local score per subject using the existing
native affine aligner. Repeated and tied occurrences are not enumerated. Scores
use declared match, ambiguity, mismatch and gap terms; they are not BLAST E-values,
expression predictions, sample identity or cloning verification. The report
retains certain/ambiguous counts, gap statistics, original subject spans and
original query coordinates after reverse complementation. Circular hits can have
two subject spans.

An exact canonical full-query substring attains the mathematical upper bound:
each query base contributes at most the positive match score and every gap cost
is negative. Such a hit proves the best score using a bounded substring scan,
without allocating a DP matrix or enumerating circular cuts. Ambiguous queries
use the native aligner and never gain invented certain matches. A circular
reference is extended by at most query length minus one, and only when the query
fits within one lap.

`best_score_proven` records score completeness for each subject. Untested origins,
cell/memory/work bounds, cancellation or incompatible alphabets keep it false
unless an upper-bound score already proves the maximum. Overall `complete` means
all requested subjects' best scores are proven under this declared profile;
it does not mean every occurrence was enumerated. Hit truncation is a separate
flag. Returned hits sort by score, then observed certain/paired counts and stable
source key. Tied alignments' alternative coordinates remain outside this profile.

## Bounds and current evidence

Inputs contain at most 128 subjects and 8 MiB total sequence bases, with aggregate
native model/annotation bounds checked before validation, hashing or rotation.
Options permit 1–64 circular origins, at most 64 million cells per alignment and
128 million total work units. Defaults use 8 million cells per alignment and
64 million total units. Work admission counts scan input lengths and complete DP
cell counts; it is not a wall-clock estimate. Cancellation is checked before each
scan and DP invocation. Rejected matrix admission does not allocate the matrix.

A cell-count multiplication that exceeds WASM32's integer range produces the
same incomplete subject result as the native cell cap. Oversized input collections
are rejected rather than silently shortened. Results never alter sources,
annotations, history, candidate plans or verification status.

Independent tiny affine-score fixtures enumerate every origin on small circular
subjects and both query orientations. Other cases cover literal origin-spanning
and reverse hits, indel/ambiguity statistics, one-lap refusal, stable ranking,
source identity, RNA, incompatible alphabets, cancellation and resource limits.
These fixtures establish the bounded core profile, not end-to-end application
or whole-genome search acceptance.

## Required application integration

The worker boundary must bind the current report version/context and bounded
subject/result shape; it must not accept a producer's completeness claim as a
cloning verification result. Web controls must capture exact query/subject event
bindings plus the workspace epoch, so equal-content source changes still cancel
stale results. Source navigation uses the retained original spans and never edits
sequence origins. MCP must resolve immutable document references and return
bounded analysis evidence, with no candidate/material installation implied.
