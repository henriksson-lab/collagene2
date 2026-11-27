# Antibody numbering implementation

`collagene-antibody-numbering` is a new independent MIT Rust/WASM module for exact
position labels, the published IMGT CDR3 allocation rule and an embedded profile-
alignment kernel. These are numbering components, not automatic antibody or
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


## Embedded human V-gene references

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
Runtime loading performs no network or executable calls. Corpus version 1 rejects
unsupported internal versions; its content revision binds the complete retained
reference data and source/exclusion metadata. Native structural validation checks
columns, identities, source ownership and region agreement. A separate integration
fixture checks every generated protein against the existing Rust translation
engine, independently of the Python generator's codon dictionary.

These are germline V-prefix facts. Unlabeled CDR3 tails are excluded; there is no
rearranged CDR3, J/FR4 assignment, nonhuman coverage, Kabat/Chothia correspondence
or calibrated classification model in this corpus. Reference metadata is distinct
from an automatic query-domain assignment. Automatic detection still requires
profile construction/calibration, ambiguity and alternative-domain handling,
numbering reconstruction and independent biological validation, followed by the
worker/MCP/protein/alignment tracks. The corpus does not complete those gates.

[Source attribution, license and rebuild instructions](../crates/antibody-numbering/data/ogrdb-human-v/ATTRIBUTION.md)
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

`v_mapping::map_query` returns all three locus alignments independently. It does
not select a chain or species by comparing scores from differently sized profiles.
Above the declared score filter, a unique path through a fully resolved query can
produce candidate reference correspondences, deleted model positions and query
insertion indices. Inserted residues receive no invented numbering labels. Any
ambiguous residue anywhere in the source conservatively withholds the projection;
tied optimal paths also withhold it. Alignment bounds and paths remain available
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
J/FR4 and rearranged CDR3, accepted alternative schemes and nonhuman scope,
alternative-domain handling, region derivation and application tracks remain open.
