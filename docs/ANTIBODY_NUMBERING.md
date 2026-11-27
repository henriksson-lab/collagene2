# Antibody numbering implementation

`collagene-antibody-numbering` is a new independent MIT Rust/WASM module for exact
position labels, published IMGT CDR3 allocation, profile alignment and conditional
V/J candidate reconstruction. These are numbering components, not automatic antibody or
variable-domain detection. No worker, UI, MCP integration, full-domain numbering,
Kabat/Chothia assignment or automatic CDR/framework derivation is implemented yet. Those remain part of F34.

## Labels and conditional CDR3 allocation

A `Position` retains an integer number and an explicit insertion kind: none,
uppercase letters or a positive decimal insertion index. Labels such as `111.1`
and `111.10` remain distinct. Decimal suffixes are indices, never floating-point
fractions. Syntax parsing is strict and cannot certify a numbering scheme or a
biological position. It rejects leading-zero aliases, lowercase insertion letters,
zero decimal indices and malformed labels. No generic ordering is supplied:
source order must follow the actual scheme, not lexical or floating-point sorting.

`imgt_cdr3_positions` follows the published [IMGT CDR3 correspondence chart](https://www.imgt.org/IMGTScientificChart/Numbering/CDR3-IMGTgaps.html).
For shorter loops it leaves the specified central gaps. Longer loops receive
alternating extra positions on the 112 and 111 sides; source order traverses the
111 insertions upward and the 112 insertions downward. The supported length range
is 2–2,048 residues as a resource-bounded rule allocator, not a claim that every
such loop is biologically realizable. This function does not renumber other loops
or translate labels between schemes.

`number_anchored_imgt_cdr3` takes an exact immutable protein source and explicitly
provided zero-based Cys/J-anchor coordinates with provenance. It checks literal C
and F/W-G-X-G anchor context, source bounds, loop stops and length, then allocates
labels to the intervening source residues. The flanking anchors are excluded from
the CDR3 residue list. Ambiguous conserved anchor identities refuse. This is
conditional on the user's anchor interpretation: those letters alone cannot prove
that the positions are actual IMGT 104/118 or that an antibody/domain exists.

Source, anchor assumptions, rule data and results remain bound together. The
streaming source hash matches the shared molecule revision; source annotation
changes also invalidate prior contexts. `ConditionalCdr3::validate_for` reconstructs
from independently supplied source and anchors. Current engine is
`antibody-numbering-1`, with no older internal reader. The conditional source
profile is capped at 2,000 residues and uses the shared model/text/location bounds.

## Full-domain implementation still required

Numbering, domain classification and CDR definitions are separate layers. An
automatic detector must derive candidate variable-domain alignments against
explicit independently checked references, preserve ambiguous/partial outcomes,
and assign actual labels before deriving regions. It must support the accepted
scheme/chain/species scope, independently verify assignments and expose the
numbering/region tracks on proteins and alignments. Raw residue offsets and motif
matches cannot replace these requirements.

The [IMGT chart](https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html)
also warns that some cross-numbering correspondences cannot be converted
automatically. CDR1/2 positioning has additional context-specific rules; the
CDR3 allocator must not be reused as an unverified general loop-numbering rule.
Future reference manifests must retain source bytes/hashes, explicit numbering
scheme evidence, authorship and dataset terms. Depositor residue numbers alone
will not be assumed to be Kabat, Chothia or IMGT numbering.

The [wwPDB archive policy](https://www.rcsb.org/pages/policies) provides a CC0
source for future structure-derived fixtures. That policy applies to archive
data, not every integrated external annotation. No IMGT database dump, third-party code or trained model is bundled. The human
reference corpus below uses separately attributed CC0 OGRDB data; conditional
anchor-rule fixtures remain entirely authored examples.

Four native fixture groups and the optional-schema WASM check pass. Fixtures
include literal published loop allocations, distinct long insertion labels,
a 2–2,048 length/uniqueness/boundary sweep, source/provenance/report tampering,
resolved and unresolved anchors, interrupted loops and resource/coordinate refusal.
These are not evidence of automatic detection, full-domain numbering or integrated
GUI/platform acceptance.

## Profile-alignment kernel

`profile::align` performs bounded semiglobal affine alignment against a declared
integer-scoring profile. Source flanks are free; the whole profile must be consumed.
Match, deletion and insertion states have explicit transitions, with no direct
insertion/deletion edges. Required match states cannot be deleted, and insertion
permission is specified after each model state. Stops block emitting states; a
new alignment may start after a stop. The profile must have at least one required
state, avoiding empty all-deletion matches.

The native profile version is 1 and kernel engine is
`antibody-profile-alignment-1`. Emission columns use the shared 22-residue ordering
in `collagene-model`. Null forbids an emission; finite weights and nonpositive
affine gap costs are bounded integers. No unvalidated third-party trained profile
is bundled. A profile name/provenance string does not certify its biological
calibration or provide a numbering scheme.

B/Z/J/X alternatives come from the shared molecular model, also used by liability
screening. The engine computes a conservative lower bound on the best score using
minimum emission scores, and an upper bound using maximum scores. The lower bound
need not be tight; absent lower paths can coexist with possible upper paths. The
upper bound is attainable across residue assignments. Threshold status describes
only these score bounds under the supplied model, not antibody identity or a
numbering certificate. Neither probabilities nor biological risk scores are emitted.

Each bound retains a deterministic representative state path with original source
and model indices, plus an explicit flag for multiple optimal paths. A tied
representative must not be silently promoted to a unique numbering map. Endpoint
ties choose the earliest source end; state ties choose match, deletion, insertion
order. Insertions and deletions remain in the path instead of disappearing into
an offset. Full automatic numbering still needs validated profiles, classification,
path-to-label rules, ambiguity handling and independent biological fixtures.

Profiles are limited to 512 match states, sources to 2,000 positions and the
caller-selected grid budget to two million cells. Every cell holds three states;
the same bounded matrix is reused for both bounds. Emissions are cached by input
symbol. Cancellation produces an error rather than a partial successful report.
Source, complete profile, algorithm/topology interpretation and limits bind the
report; `profile::Report::validate_for` reconstructs all scores, paths and tie flags.
It validates computation under the profile, not the model's biological accuracy.

Five additional native fixture groups check explicit affine paths, required and
optional states, tied endpoints, stop barriers, thresholds, impossible emissions,
ambiguity bounds, current-version/budget checks, cancellation and report tampering.
An independent recursive path-language oracle covers all 340 canonical inputs of
lengths 1–4 for three small profiles (1,020 comparisons); concrete residue expansion
checks ambiguous score bounds. Together with the four numbering and seven affected
liability fixtures, 16 native cases pass. Optional-schema WASM checks cover both
consumers of the shared amino-acid alternatives.


## Embedded human V/J references

The `references` module loads a version-pinned CC0 corpus from OGRDB releases
G00060 v10, G00061 v5 and G00062 v4. It retains 354 provider-annotated functional
human V alleles: 194 heavy, 74 kappa and 86 lambda. Explicit provider IMGT labels
cover positions 1–104. Aligned nucleotides, derived protein columns, allele/source
identity, CDR1/2 boundaries and complete inclusion/exclusion accounting travel
together. Gaps are preserved; six alleles with atypical conserved cysteines remain
present and flagged. They must not be silently corrected or admitted as canonical
conserved-anchor references.

The source API responses are pinned by SHA-256. An authored generator verifies
source identity, frame, codon gaps, prefix agreement and label/CDR agreement.
Runtime loading performs no network or executable calls. Corpus version 2 rejects
unsupported internal versions; its content revision binds the complete retained
reference data and source/exclusion metadata. Native structural validation checks
columns, identities, source ownership and region agreement. A separate integration
fixture checks every generated protein against the existing Rust translation
engine, independently of the Python generator's codon dictionary.

These are germline reference facts. Unlabeled V CDR3 tails are excluded; there is no
rearranged CDR3, query-domain assignment, nonhuman coverage, Kabat/Chothia correspondence
or calibrated classification model in this corpus. Reference metadata is distinct
from an automatic query-domain assignment. Automatic detection still requires
profile construction/calibration, ambiguity and alternative-domain handling,
numbering reconstruction and independent biological validation, followed by the
worker/MCP/protein/alignment tracks. The corpus does not complete those gates.

[Source attribution, license and rebuild instructions](../crates/antibody-numbering/data/ogrdb-human-vj/ATTRIBUTION.md)
· [Provider content license](https://ogrdb.airr-community.org/render_page/licensing_statement.html)
· [AIRR numbering/delineation contract](https://docs.airr-community.org/en/latest/datarep/germline.html).

## Query-to-reference V-prefix mapping

`v_mapping::models` constructs three explicit human reference models (heavy,
kappa, lambda) from the current corpus. The six atypical Cys23/Cys104 alleles are
excluded from this canonical-anchor modeling policy, with their identities
retained. The other 348 alleles yield 140, 60 and 72 distinct aligned protein
sequences, respectively. Identical protein sequences count once, irrespective of
how many synonymous alleles encode them. This is reference-sequence weighting,
not an estimate of repertoire abundance or a biological prior.

Each populated column retains its actual provider position. All-gap columns
have no emitting model state and stay explicit as omitted reference positions;
subsequent positions are not renumbered. Per-column counts include gaps and every
residue. For an observed residue, the integer score is `observed_reward +
floor(frequency_bonus * count / occupied_sequences)`; unobserved residues receive
the separately declared penalty. Required Cys23/Cys104 states emit only C. Other
populated states permit deletions and insertions according to the affine kernel.
These are authored heuristic scores, not trained likelihoods or calibrated
antibody evidence. Callers supply all weights and their score threshold explicitly.

`v_mapping::map_query` searches all three loci independently at every full-profile
source endpoint. It does not select a chain or species by comparing scores from
differently sized profiles. Above the declared score filter, a unique path at an
endpoint through a fully resolved query can
produce candidate reference correspondences, deleted model positions and query
insertion indices. Inserted residues receive no invented numbering labels. Any
ambiguous residue anywhere in the source conservatively withholds the projection;
tied optimal paths at that endpoint also withhold it. Distinct endpoints remain
separate candidates. Alignment bounds and paths remain available
for inspection. Even a unique, high-scoring correspondence does not certify an
antibody, complete domain, actual IMGT numbering or a CDR/framework assignment.
The current model is not an independently calibrated domain detector.

Model construction, corpus identity, explicit scoring, profile-engine semantics,
source revision and limits bind the mapping report. Model validation reconstructs
from the current embedded corpus; report validation reconstructs every alignment
and projection. The cell budget covers all three model grids together, not three
separate uses of the caller's budget. Cancellation never returns a partial
successful report. Source and per-model bounds remain those of the native kernel.

Native fixtures check independent literal reference counts and missing columns,
real V-prefix source anchors, all-locus evidence, duplicate best matches,
ambiguity, declared thresholds, budgets, cancellation and tampering. This does not
replace independent biological calibration and full-domain scheme validation.
Validated full-domain assignment, accepted alternative schemes and nonhuman scope,
alternative-domain handling, region derivation and application tracks remain open.


## J/FR4 reference facts

The current corpus also retains all 23 provider-annotated functional human J
alleles from the pinned releases (7 IGH, 7 IGK, 9 IGL). Original sequence, coding
interval, coding DNA and raw provider frame/anchor values survive the transform.
For these OGRDB responses, the conserved anchor is verified as a coding-sequence
offset and its phase agrees with the complete-codon start. This is an explicitly
checked provider-specific convention, not a generic AIRR import rule.

Complete codons starting at the conserved W/F yield FR4 reference proteins. The
canonical motif is checked against the original DNA translation. Labels start at
published IMGT position 118 and run sequentially across complete codons: through
128 for these heavy references, 127 for these light references. These J labels
are authored rule-derived facts, not provider-supplied per-column labels. The
terminal one-nucleotide partial codon remains separate and untranslated; it is
not assigned an FR4 amino acid or combined with invented constant-region bases.
[IMGT FR4 definition](https://www.imgt.org/IMGTeducation/Tutorials/3Dstructure/_UK/DOMAINs/fleches.html)
· [Coding/splice boundary](https://imgt.org/IMGTScientificChart/Numbering/IMGTnumberingCR.html).

The single current corpus representation is version 2 / engine
`antibody-reference-corpus-2`. Old internal versions are rejected. V-reference
models and tests consume the updated corpus directly; there is no V-only adapter
or fallback. Counts and source-exclusion ledgers cover V and J separately. Native
fixtures check every J translation with the existing Rust engine, including the
retained untranslated nucleotide, and refuse invalid raw coordinates, anchor
phases, label extents and terminal metadata.

These are J reference facts. The conditional pairing section below connects them
to query V paths. Full scheme validation/calibration, alternative-domain handling
and worker/MCP/application tracks remain required.


## Bounded endpoint search

`profile::align_endpoints` returns the best full-profile alignment **at each source
end** that meets the declared attainable-score threshold. It uses the same matrix
fill, transitions and traceback as the general single-best query. A single matrix
is reused for upper and lower bounds; traceback retention does not require one
matrix per endpoint. Candidate endpoints are retained in ascending source order.
Reachable and score-qualified endpoint counts cover the whole input even after
the retention cap is reached. A zero cap requests counts only. A truncated result
must not be interpreted as absence of additional qualifying endpoints.

The endpoint engine is `antibody-profile-endpoints-2`, with at most 128 retained
endpoints per profile. Each endpoint has its own score bounds, source/model path
and tied-path flag. Repeated independent sequence blocks therefore retain separate
endpoints. Multiple optimal starts or state paths **at the same end** remain tied;
they are not exhaustively enumerated. This mode also does not enumerate suboptimal
starts at a shared endpoint. Its counts describe score-qualified endpoints, not
biologically confirmed domains or all possible numbering assignments.

V mapping now uses engine `antibody-v-reference-mapping-4`. Its former single-best
representation is replaced directly: one endpoint report per locus and explicitly
end-bound projections. The caller supplies an aggregate grid-cell budget and a
per-locus retention cap. Counts, truncation, thresholds and exact model/source
context survive report validation. Unsupported mapping engines reject without
an adapter or fallback. Biological calibration, full-domain numbering and region
derivation, accepted scheme/species coverage and worker/MCP/application tracks
remain open.

Native endpoint fixtures compare all 1,020 small input/profile pairs at every
source end against independent complete path enumeration, including terminal
insertions and deletion paths. They also check concrete ambiguity bounds,
thresholded and counts-only results, retention limits, reconstruction and
cancellation. Tandem real-reference V-prefix fixtures retain both endpoints and
exercise truncation without turning duplicate blocks into an unreported global tie.


## Conditional V/J candidates and CDR3 reconstruction

`domain_candidates::find_candidates` connects retained V endpoints with an
invariant terminal C104 anchor to downstream complete FR4 windows in the same human locus profile.
Callers explicitly declare V scoring, CDR3 length limits and an FR4 substitution
budget (0–3). No J indels or partial FR4 tails are accepted by this finite profile.
The locus-specific W/F anchor and G119/G121 motif positions remain literal, even
when substitutions elsewhere are permitted. Stop-containing loops/tails refuse.

Identical FR4 protein references are grouped with every original allele ID.
Different allowed reference tails at one source interval remain explicit J-match
alternatives with source-indexed substitutions. Different downstream J starts
remain separate candidate intervals; the algorithm does not select a closest J
sequence or one preferred domain. Every retained candidate receives conditional
IMGT CDR3 labels from the existing published-rule allocator and sequential FR4
labels from the actual query residues. V correspondence gaps/unassigned residues
remain in the linked V evidence; they do not acquire invented labels. Repeated
hashing of the full immutable source is avoided during loop allocation.

Search counts are scoped to retained endpoints with an invariant C104 anchor.
Ambiguous sources withhold that anchor and are reported as blocked. Tied internal
V positions remain unprojected, but their fixed terminal anchor permits conditional
CDR3 pairing; each candidate retains its V projection issue. V endpoints without a
qualifying complete J window remain listed. The output independently reports V
retention truncation and pair retention truncation. A zero pair-retention cap still
counts all qualifying intervals in the retained-anchor scope. A zero V cap searches no
V/J pairs and retains the V truncation flag, rather than claiming global absence.
The algorithm does not enumerate suboptimal V paths or infer omitted endpoints.

At most 64 candidate payloads and two million complete-FR4 scan windows are
allowed; callers can choose smaller limits. Exceeding work bounds or cancellation
returns an error, never a partial successful report. Native engine
`antibody-vj-candidates-4` binds the exact source, current V/reference/endpoint and
CDR3-rule engines, options and limits. Report reconstruction checks all candidates,
coordinates, loop labels, FR4 residues, J alternatives, counts and truncation flags.

Fixtures use authored joins of attributed CC0 reference prefixes/tails, not
naturally observed biological rearrangements. They cover heavy/kappa/lambda
profiles, tandem constructs, identical-tail allele groups, competing J sites,
substitutions, scope/stop/partial-tail refusal, bounded/counts-only retention,
source tampering and cancellation. They validate conditional reconstruction,
not biological detector accuracy. Independent biological calibration, complete
scheme/region assignment and application integration remain required.


## Invariant anchors and external assessment

A required, non-insertable final V model state represents C104 and emits only C.
At a fixed source endpoint every state path therefore assigns the same source
residue to this terminal reference position, even if internal paths tie. The
current mapper validates this topology and retains the C104 match among
`invariant_correspondences`, independently of its full projection. Tied V paths still receive no invented unique V-position
map. Conditional V/J pairing can use the invariant anchor while retaining
`v_projection_issue`; its `representative_v_start` explicitly describes the stored
representative rather than proving a unique domain start.

The development assessment panel contains 30 deposited whole protein chains from
13 included RCSB entries, with 12 antibody-chain and 18 non-antibody identities fixed before
algorithm inspection. It excludes four explicitly recorded short/ambiguous
entities. Antibody controls include real human/humanized Fabs and intact IgG;
non-antibody controls include TCR, MHC, CD4/CD8, beta-2 microglobulin and other
proteins. Only archival sequence/molecule facts are retained; third-party numbering
annotations are not copied. This is a chain-identity assessment, not independent
residue-level numbering truth or population-wide calibration.

The initial run found four antibody chains without candidates because internal
V ties unnecessarily blocked the fixed C104 anchor. The current anchor evidence
separates that ambiguity from conditional loop reconstruction. Native regression
fixtures retain the V tie and check the independently sourced sequence case,
source anchor binding and reconstructed CDR3. Broader biological calibration and
complete scheme/region validation remain required.

[Panel attribution, source manifest and reproduction](../crates/antibody-numbering/data/rcsb-assessment/ATTRIBUTION.md).


The assessed options retain V weights 5/8/−6, affine gaps −12/−2, minimum V score
300 and CDR3 lengths 2–60. Only the declared FR4 substitution allowance varies.
Observed with candidate engine 2 / V-mapping engine 4 on this fixed panel:

| Allowed FR4 substitutions | Antibody chains with a matching chain-class candidate | Non-antibody chains with candidates |
| --- | --- | --- |
| 0 | 10 of 12 | 0 of 18 |
| 1 | 12 of 12 | 0 of 18 |
| 3 | 12 of 12 | 0 of 18 |

The exact-tail profile misses both 1HZH chains; their observed FR4 sequences differ
from the reference tails. No V or pair retention cap is reached in these runs.
The current query-loop batch passes all 14 selected native cases and focused
optional-schema WASM checks,
including source-bound invariant-anchor and observed-sequence regressions. These
results do not validate lambda detection, nonhuman coverage, arbitrary mutation
loads or residue-level numbering. The panel is small, deliberately selected and
contains related sequences; it is not a population performance estimate.


## Invariants over every optimal path

The endpoint kernel now stores all tied predecessor states in the dynamic-
programming graph. A bounded backward traversal visits the optimal-state union;
a model position is invariant only when every optimal path matches it to the same
source index and no optimal path deletes it. The deterministic representative and
its tie flag remain available separately. This computes exact intersections of
path correspondences, not a chosen path with an ambiguity label removed.

Invariant summaries are withheld for residue-ambiguous sources: the graph of an
upper score bound alone cannot prove correspondence over every possible residue
assignment. The independent `max_trace_states` budget limits graph states across
all retained endpoints and, in V mapping, across all three loci. At most two million
states may be visited. Budget exhaustion or cancellation returns an error; a
counts-only request needs no invariant trace work. The report retains exact work
counts and binds the budget into its context. Current endpoint/V/candidate engines
are 2/4/4; callers and fixtures use the new representation directly.

An independent complete-path oracle now checks both endpoint scores and match
intersections for 1,700 small input/profile configurations, including zero-cost
gaps, terminal insertions, deletions and multiple starts. Existing single-best
score/path checks remain unchanged. The 30-chain assessment retains its previously
reported outcomes under the same scoring and substitution options.

## CDR1/2 central-gap rule allocators

`imgt_loops::positions` implements the published V-domain display convention:
CDR1 occupies 27–38 with lengths 5–12; CDR2 occupies 56–65 with lengths 0–10.
Short loops keep the extra residue on the left when their length is odd, leaving
central gaps. It rejects unsupported lengths and invents no insertion labels.
This differs from rearranged CDR3 allocation and from moving/storing germline
positions for other analyses. [IMGT V-domain chart](https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html).

Native literal chart fixtures cover all supported lengths; both loop placements
also agree with all 354 provider-labeled V references. The rule engine is
`imgt-cdr12-central-gaps-1`, with an explicit provenance revision. Query reconstruction below applies these allocators only when the corresponding
framework anchors are independently established. Full framework numbering,
alternative schemes, broader biology validation and application tracks remain
unfinished.


## Source-bound query CDR1/2 reconstruction

`imgt_loops::number_anchored_loop` validates an explicit immutable protein source,
provided framework boundary indices and provenance, then reconstructs a conditional
CDR1 or CDR2 report. Anchor residues are excluded. These boundary positions do not
have fixed amino-acid identities, so coordinate validation does not independently
prove their biological interpretation. Stops at either boundary or inside the
loop refuse; unsupported lengths refuse. A zero-length CDR2 has an explicit empty
residue list and equal start/end, distinct from an unresolved loop. Source, kind,
anchors and published-rule revision bind report reconstruction.

Candidate engine 4 now carries CDR1 and CDR2 outcomes alongside CDR3/FR4. It uses
only reference positions 26/39 or 55/66 present in every optimal V path. Complete
boundary evidence yields source-bound loop labels in the published central-gap
layout; missing invariant boundaries or unsupported lengths yield typed unresolved
outcomes. The existing CDR3 candidate and its other evidence remain available.
No representative-path boundary is silently substituted. Shared internal
allocation reuses the already validated source revision instead of hashing the
whole molecule again for each candidate loop.

The assessment example includes both loop outcomes for every candidate. Native
fixtures check literal query coordinates/residues/labels, empty CDR2, boundary and
stop refusal, provenance/source binding, report tampering and an engineered long
loop whose unresolved status does not erase CDR3. Full framework numbering,
Kabat/Chothia and other accepted definitions, broader validation and application
tracks remain required.


On the fixed external panel with one or three allowed FR4 substitutions, the
12 matching-chain candidates contain 23 reconstructed CDR1/2 outcomes and one
unresolved outcome: 1HZH entity2 CDR1 lacks an invariant position39 boundary.
That case retains its conditional CDR3 and explicit diagnostic. Candidate presence
and negative-control results remain unchanged. These conditional loop results do
not establish independent residue-level numbering accuracy across the panel.

## Conditional source-coordinate tracks

`domain_tracks::analyze` computes V/J candidate evidence once and derives one
track for every retained candidate. The report carries that exact evidence,
source revision, published-rule revision and analysis context. Numbering scheme
and region definition are separate typed choices; currently only IMGT/IMGT is
implemented. Unsupported choices are rejected rather than converted by offsets.

Framework assignments come only from reference correspondences shared by every
optimal V path. CDR1/2 storage labels are excluded: reconstructed loop rules supply
their display labels, and an unresolved loop receives no replacement labels.
Conditional CDR3 and matched J/FR4 complete the supported assignments. Each
assignment retains its framework-reference, conditional-loop or J-reference
origin. Regions use the published boundaries FR1 1–26, CDR1 27–38, FR2 39–55,
CDR2 56–65, FR3 66–104, CDR3 105–117 and FR4 118–128.
[IMGT V-domain chart](https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html).

Every source residue from the earliest supported framework correspondence
through the candidate FR4 end has an absolute source-coordinate row. Unassigned
rows remain explicit and counted. `invariant_position_one` is present only when
reference position 1 has an invariant match; its absence leaves the N-terminal
boundary unresolved. Earlier source residues are not asserted to lie outside the
domain, and the representative alignment start is never used as a substitute.
Empty CDR2 remains explicit in the attached candidate evidence. Decimal CDR3
labels remain exact and in source order, including decreasing 112 insertion
suffixes; there is no floating-point sorting or invented framework insertion.

Tracks are conditional reference/rule assignments, not a calibrated classifier or
independent biological numbering certification. Candidate truncation and unresolved
loop diagnostics remain in the same report. Source/options/report reconstruction,
source bounds, cancellation and existing candidate caps apply. Protein/alignment display beyond the inspector, alternate schemes/definitions and
broader independent numbering validation remain required.


## Worker, Biology inspector and MCP

The current track engine is `antibody-domain-tracks-2`. Its request context binds
source, corpus, scoring/loop/FR4 options, work/retention limits and current engines
and rules without running alignment. `validate_binding` checks current context and
structural track/evidence consistency for trusted local worker replies: source
coordinates, labels, regions, assignment origins, repeated overlapping source
letters, candidate identity and retained-count flags. It does not establish DP
scores or authenticate source letters. Imported or independently supplied reports
still require `validate_for` against an independently resolved actual source.
There is no decoder for the superseded internal track context.

Worker schema 118 adds job 195/result 185 with ordinary message bounds. The Biology
inspector runs analysis in a cancellable worker, offers explicit analysis settings,
retained-candidate selection, loop diagnostics, paged residue rows, exact source
selection and full JSON download. Document ID, event, revision, workspace epoch
and options/limits pin every request and navigation action. Equal-content foreign
documents cannot retarget output. Source/settings changes invalidate old results.

Native MCP offers `list_antibody_profiles` and `analyze_antibody` against immutable
protein documents, explicit options and limits. At most 20 candidate summaries are
returned inline; a full report artifact preserves all retained evidence. These
read-only tools do not install documents or materials. The inspector and MCP expose
the conditional human-reference scope directly. They do not complete alternate
schemes, independent biological validation, alignment-row projection or full
whole-product acceptance.
