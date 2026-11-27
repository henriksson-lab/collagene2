# Sequence similarity across selected documents

The native `sequence_similarity` module, binary worker, Compare panel and
`search_sequence_library` MCP tool share one bounded collection-search profile.
The scientific engine is `sequence-similarity-2`; the current worker schema is
115 (job191, result181). Previous worker versions are rejected directly.

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

## Application workflow

Open **Compare → Search similar sequences in this workspace**. The active
DNA/RNA document is initially pinned as the query. **Use active document as query**
changes it explicitly; navigating a match does not retarget the pinned query.
Choose the document scope, review **Scoring and search limits**, then search.
Matching documents are selected by default, excluding the query itself; selecting
the query explicitly is valid. Missing documents and more than128 selections
require scope review rather than silent truncation.

Results show score, paired-query coverage, strand and circular seams. Click a
match to navigate its exact source event and original coordinates. Results page
20 rows at a time. **Search evidence and limits** retains unsearched/incomplete
rows; **Download search evidence** includes report and exact event captures.
Changing sources, event identities, options, scope or workspace epoch invalidates
results and terminates pending work. Cancel/close also terminate the isolated
worker. Biological documents and origins are never edited by navigation.

The borrowed native preflight runs before selected molecules are cloned. Worker
responses check the current engine, source/options context, bounded dimensions,
statistics, ranking and completeness consistency. Web captures additionally bind
every selected event, alphabet, topology and length; equal-content events cannot
retarget a reply. Completeness is evidence about this declared search profile,
never molecular verification authority.

MCP resolves immutable `DocumentRef` query/subject handles through its session
store and applies the same borrowed preflight before cloning. The tool returns
bounded summaries plus a complete JSON artifact, with explicit incomplete rows.
Stale handles fail before analysis. No document, candidate, material or cloning
run is installed. This profile searches explicitly loaded workspace/session
subjects; a server-wide or genomic indexed similarity service is not claimed.
