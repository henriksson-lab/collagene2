# Junction edits and coherent overlap-primer proposals

`assembly_junction_edit` separates two actions: freeze an independently declared intended reference, then propose PCR-overlap primers that could prepare it. Neither action copies a simulated assembly answer into its own expected target. Original source molecules, source revisions and source coordinates remain the inputs to every generated PCR operation.

An `Edit` identifies junction `i` after row `i`, with `remove_left` bases removed from that upstream oriented core, `remove_right` bases removed from the next oriented core, and literal `inserted` DNA between them. Linear assemblies have `n−1` junctions; circles additionally have the last-to-first junction. Each junction may occur once. Incoming and outgoing trims together must retain at least one original base in every row.

For a forward source contribution, the incoming trim increases its original start and the outgoing trim decreases its end. For a reversed contribution, the outgoing trim increases the original start and the incoming trim decreases its end. Circular source coordinates remain one-turn intervals; a trimmed start passing the original source origin is normalized together with its unrolled end. Selection orientation is applied after extraction, not by editing a reverse-complemented source document.

`build_intent` first builds the existing pure assembly intent from the adjusted source cores. It then inserts only the declared literal DNA at the original concatenation boundaries, in descending order, using native sequence editing and annotation review rules. The resulting reference is the concatenation of each adjusted oriented core and its declared outgoing synthetic suffix. Circle-closing inserts are appended after the final core, preserving the first adjusted core as the intended origin. No PCR, ligation or overlap simulator is called.

The returned `Intent` retains the original base request, exact edits, adjusted core-only request, source revision pins, target and source contributions. A contribution's output range includes its core **plus its owned synthetic suffix**, so output ranges cover the reference without gaps. Its source start/end still identify only the original source core. Synthetic insert annotations explicitly mark the added bases, require review and make no inherited-function or synthesis claim. Existing source annotations are projected by the shared engine; clipped or affected annotations retain review warnings, while complete original source evidence stays in the original input documents. Targets carry a restrictive reference-only molecular-state marker.

`design` uses the adjusted original-source intervals with the existing in-house anchored primer designer. At each outgoing junction, the intended DNA added to the upstream oriented core is:

```
upstream literal insert + next adjusted oriented core prefix(overlap)
```

The reverse complement of that complete addition becomes the upstream original reverse primer's 5′ tail for a forward row, or its original forward primer's 5′ tail for a reversed row. The other tail remains empty. The last linear row adds no tail. Circular closure extends toward the first adjusted prefix. Each row is analyzed and ranked with its actual complete tail already present: structure screens and pair comparisons include introduced bases, while Tm retains the existing annealing-domain model and explicit calculation conditions. Names are retained as role-forward/role-reverse. User-supplied option tails must be empty; they are never silently replaced.

The resulting recipe contains existing `ExactPcr` preparations against the immutable originals and explicit ordered overlaps. The downstream recipe compiler/verifier must freshly reconstruct all preparations against `Intent.target`, retaining competing PCR geometries and uncertainty. Neither a design score, a local frame display nor a matching intended reference establishes a successful experiment.

Borrowed `intent_preflight` guards sources, edits, annotation work and output size before copying rich models. `preflight` additionally guards the existing primer-source/filter rules and charges the actual extended insert-plus-homology tails before design. Limits are 32 rows/edits; at most 100 normalized bases per insert and 200 raw ASCII bytes; at most 4096 inserted bases overall; an 8 MiB combined original-source/intended-reference DNA budget; and the existing 200-base complete ordered-oligo profile. Native per-row pair-work and aggregate 100-million-unit design/PCR scan limits remain explicit. Intent permits valid IUPAC notation; concrete oligo proposals with ambiguous cores or inserts are withheld. A long or ambiguous intent may therefore exist even when its requested preparation cannot be proposed.

Context revisions bind engine domain, actual input revision pins and the original declarations/options. Readonly `Intent::validate_binding` and `Report::validate_binding` check original/adjusted coordinate correspondence, synthetic contribution lengths, sequence receipts and complete ranked-pair/recipe/metric identity. They do not rerun primer design or restore scientific authority from imported evidence. The old pure-concatenation report validator must not be used to reinterpret the edited target.

Five focused fixtures independently specify linear and circular literal targets, reversed-source trimming, source-origin crossing after trims, annotation review/synthetic bases, original immutability, full tail metrics, actual fresh recipe reconstruction and wrong-target rejection, plus stale contexts, duplicate edits, over-trimming, ambiguity, invalid option tails and bounded-work refusal. No scientific operation, kernel, candidate schema, plan schema or verification engine version changes.

## Connected draft transaction

The desktop web assembly draft retains bounded junction edits separately from its
original pinned source slots. The central editor previews complete named oligos,
actual synthetic/homology tails, metrics and the independent intended DNA. Applying
explicitly replaces every proposed row with original-source exact PCR preparation
and updates edits/overlaps/frozen target together. Plain unedited-region suggestions
cannot overwrite active edited primers. The original region-only path remains
unchanged when no edits are declared.

Undo retains the before/applied snapshots. It is available only while declarations
and methods still match the applied snapshot, ignoring active-slot navigation; later
source/name/primer changes are never overwritten. Undo preserves the current slot
where possible and clears computed evidence. Reloaded drafts recompute intent and
analysis. Selective preservation of unaffected/mixed preparations remains pending.

Five core fixtures, two binary protocol cases, two real MCP SDK groups and seven
native staging cases passed (16 total, including four existing staging regressions).
Existing junction-inspection fixtures remain separate and unchanged. The full
130-tool stdio smoke and final web/worker WASM compiler check also passed. Browser
checks remain deferred; source/native checks do not establish GUI acceptance.
