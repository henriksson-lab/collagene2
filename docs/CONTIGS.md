# Reference-free read contigs

`core::contig` assembles immutable linear DNA read documents without a reference
sequence, an external executable, or a multiple-alignment engine. It is an
independent, conservative overlap/layout/consensus implementation, not a CAP3
port or a claim of parity with an established assembler. Conceptual background
includes [Huang and Madan, CAP3 (1999)](https://doi.org/10.1101/gr.9.9.868) and
[Rausch et al. (2009)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2732307/).
Read assembly and evidence-derived consensus are distinct from a general MSA.

## Inputs and evidence

Each input has a unique, stable role ID and an exact `Molecule` snapshot. ABIF,
SCF and FASTQ attachments retain their native score and original-call mappings.
The existing trace/quality readers validate the attachment against the current
document. Stale or malformed quality evidence does not silently become an
unscored read. Plain sequences require `allow_unscored`; the report explicitly
discloses this policy. Optional quality-window trimming and individual raw-score
thresholds are retained in the options. See [TRACE_FORMATS.md](TRACE_FORMATS.md)
and [TRACE_VALIDATION.md](TRACE_VALIDATION.md) for native score semantics.

ABIF caller scores and SCF called-base accuracies are not converted into invented
Phred probabilities. An eligible vote is a concrete A/C/G/T call satisfying the
declared policy, not a calibrated probability of correctness. Repeated native
attachments or FASTQ records, and canonical identical unscored sequences, are
disclosed as duplicate groups and contribute at most one group per column.
Distinct records do not establish distinct physical molecules.

## Declared overlap and layout model

Affine dynamic programming examines each ordered pair in the enabled read
orientations. A boundary path starts with a paired call, anchors the query at
its first call, and permits a free reference prefix. Two endpoint classes are
searched separately:

* Extension reaches the reference end before the query end.
* Containment reaches the query end on a paired call, at any reference endpoint.

For each endpoint, the engine considers affine-optimal paths. It selects the
highest-scoring threshold-qualified candidate in each class: a positive score,
enough matching **eligible** concrete calls, and the declared error fraction.
It does not enumerate every suboptimal path or every possible assembly. Error
columns include substitutions, ambiguous pairs and indels. Co-optimal path and
endpoint counts saturate at two and prevent arbitrary unique joins. When a
representative fails a filter, a bounded traversal of optimal predecessor paths
can prove that no alternate qualifies. If an uninspected or qualifying alternate
path could pass,
the affected search is explicitly incomplete. Score and eligible-anchor bounds
may prove such a path ineligible; a saturated count alone is not that proof.

Layout joins require reciprocal unique best extension choices; when both
orientations are searched, the reverse-complement counterpart must agree.
Ordinary transitive overlaps remain inspectable even when a shorter extension
determines the path. Contained and duplicate records are retained as evidence.
A contained read is recruited only when its candidate host placements agree on
the complete call-to-column mapping. Branches, tied placements, repeat borders
and cycles remain explicit alternatives or unresolved singletons. No arbitrary
circular assembly or repeat-copy count is inferred. Progressive insertion
columns can reuse compatible existing insertion intervals; incompatible complex
layouts remain unresolved instead of forcing an MSA.

`layout_unique` refers only to this bounded, declared model. It does not prove a
unique physical molecule, biological origin, or exhaustive assembly search.
`report.complete` says the declared analysis finished without preparation,
resource or unresolved path-eligibility failures. A complete report can contain
ambiguous layouts. `contig.complete` additionally requires a unique declared
layout and resolved evidence at every column.

## Consensus and coordinates

Every emitted column retains concrete source observations. Eligible groups must
agree unanimously; majority voting never hides a conflicting eligible call.
Conflicts, low-quality-only observations and insufficient support produce `N`.
A supported gap has eligible flanking calls on the same read and is omitted
from the emitted DNA while remaining in the evidence columns. Unsequenced
terminal space does not vote for a gap. Eligible base/gap disagreement produces
an explicit gap conflict. Minimum-support and optional both-strand requirements
apply to deduplicated groups. Original evidence strand is placement reversal
XOR stored source reversal: display-reversing a document cannot manufacture
opposite-strand support.

Overlap coordinates index oriented trimmed strings. Placement `call_index`
indexes the current document; `stored_call_index` indexes the original file.
Native sample indices retain both oriented and original-file coordinates.
Placement column indices refer to evidence columns, whereas
`consensus_position` indexes emitted DNA and is absent for a supported gap.
Every retained read call maps exactly once in its placement. A separate audit
reconstructs observations, raw scores, eligibility, gap flanks and consensus
from immutable source snapshots before products are returned. No reference
base, source annotation or raw attachment is copied into consensus DNA.

## Bounds, binding and reconstruction

Source preflight runs before source validation and hashing. Limits cover read
count and length, aggregate calls, rich-source bytes and structural nodes,
location depth, pair and cumulative DP/traceback work, overlap count/text,
evidence columns and output size. Exceeding a limit is an error or an explicitly
incomplete partial report; it never establishes an absent overlap or unique
layout. Current hard limits include 32 reads, 10,000 current calls per read,
100,000 aggregate calls, 16 MiB source/report size, 16 million cells per pair,
512 million cumulative work units, 4,096 edges and 100,000 columns. User options
can choose smaller work budgets. `cells_used` includes traceback work as well
as matrix cells. The optimal-predecessor matrix uses two bytes per bounded cell
(at most 32 MB at the hard pair limit), with rolling score rows. The alternate
path traversal also caps its pending stack at 100,000 entries and charges visits
to the cumulative work budget; exceeding either bound leaves the search
incomplete.

`assemble` produces source-derived reports. `batch_revision` binds exact input
roles/revisions, engine version and options after bounded preflight.
`ContigReport::validate` checks bounded historical structural consistency;
historical engine records remain historical evidence. `validate_bindings_for`
requires the current engine and exact current source/options bindings, without
repeating DP. Neither method grants scientific authority to an imported report.
`consensus_documents` performs a fresh assembly and exact report comparison
before returning source-derived molecules. Browser opening/export reruns the
worker; MCP assembles fresh before installing immutable products. Full replay
is not performed on the browser main thread.

Generated products retain engine, batch, options revision, complete options,
role/revision summary and contig identity in metadata. A [native saved assembly attachment](CONTIG_HISTORY.md) preserves the report
and exact source-event bindings as historical provenance. The product does not embed a self-referential report hash and saved
evidence does not restore current verification without fresh computation.
