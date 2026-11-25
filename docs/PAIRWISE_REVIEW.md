# Event-bound construct comparison

The comparison panel keeps its reference pinned while the main workspace opens
the query or navigates a difference. **Use active document as reference** is the
explicit retarget action. A newly chosen query uses its current event; reopening
a saved comparison instead pins both exact historical events. Equal DNA in another
document does not replace either binding.

Two explicit policies reuse the existing native affine algorithms:

- `Stored(AlignmentOptions)` compares the stored origin/orientation in global,
  local or query-in-reference mode. No strand/origin search is implied.
- `Whole(ComparisonOptions)` requests global whole-construct comparison with
  explicit query-strand, circular-origin and candidate/cell budgets. The native
  exact-identity shortcut and bounded origin search keep their existing semantics.

Coverage, compatible ambiguity, uncertain differences and unaligned flanks remain
separate evidence. Native candidate counts, requested-scope completeness, best-score
proof and limitations are displayed. A best score is not proof of a unique
alignment, equal physical topology, sample identity or experimental success.
Whole-construct policy never silently changes a local request into global mode.

## Coordinates and reachability

All alignment columns are available in 600-column pages; differences use 20-row
pages with previous/next and explicit difference-number navigation. The complete
CIGAR is paged too. Page offsets are indexed once from the actual aligned rows,
including gaps, rather than relabeling every page from source position zero.
Long difference strings have a disclosed preview; their full bases remain in
the alignment pages and complete JSON. There is no first-4,000-column or
first-200-difference-only review limit.

Text-row R/Q coordinates are alignment boundaries in the explicitly rotated
reference and oriented query. `Report::map_reference` and `map_query` convert
selected differences back to original document coordinates. An insertion's empty
reference interval is a navigable caret; a deletion similarly provides a query
caret. A reference interval crossing its circular origin retains two ordered
spans. Reverse-query mapping uses original half-open boundaries, not reversed
display offsets.

Navigation into a current one-span event selects that exact source in the main
workspace while preserving the comparison's pinned reference. Historical or
origin-spanning selections open a shared read-only reference viewer, with initial
span/caret highlights. Viewing does not rewind a head or create a document merely
to inspect old evidence. Existing source-event and workspace-epoch guards apply
before navigation, saving and JSON export.

## Native report and archive

`core::pairwise_review` wraps `alignment::align` and `compare_molecules`. It adds
content-bound options, source revisions/lengths and interpretation version
`pairwise-review-1`. Read-only `validate_binding`/`validate_sources` checks native
scope, operations/CIGAR, differences, statistics, path score, oriented source
letters and hashes without rerunning dynamic programming or establishing an
imported optimality claim. Fresh algorithm execution supplies current evidence.

Worker52 appends job119 `CompareSequences` and result111 `ReviewedComparison`.
Older version labels are rejected before decoding these payloads. The ordinary
32 MiB limit remains; no workspace-size whitelist is expanded. Existing alignment
job/result tags remain supported. Candidate11, plan9, science-14 and archive1
are unchanged.

`core::pairwise_history` stores version1 `pairwise-review-v1` workspace
attachments with exact reference/query events and the full report. Capture
validates the privately computed report and stages only the new attachment; it
does not duplicate the scientific alignment run or mutate any head. Reports and
JSON export use the existing native entry bound. Content-addressed attachment
hashes distinguish records bound to equal-DNA aliases. Invalid records are an
explicit error, not silently reassigned sources or an omitted first-N list.

Saved comparisons reopen as historical results, never with a restored fresh
flag. **Compare pinned versions** explicitly reruns the current engine on those
retained inputs. Saving current evidence requires an unchanged source/options
context; cancellation and late worker responses cannot restore an obsolete result.
Normal autosave and native workspace export retain comparison attachments.

No external alignment executable, multiple-alignment engine, read mapper or
BAM/index integration is introduced. Integrated reference/read/chromatogram review
is implemented by [read checks](READ_CHECKS.md) and [explicit review](READ_REVIEW.md).
Computed pairwise search, aligned context and portable exports remain follow-ups
([current decisions](BENCHLING_WORKFLOW_DECISIONS.md)). Browser acceptance stays
deferred under the current instruction.
