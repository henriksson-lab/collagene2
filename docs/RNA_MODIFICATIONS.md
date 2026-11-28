# Declared RNA chemistry

`collagene-rna-chemistry` provides an independent MIT Rust analysis core for
source-bound RNA base, sugar, backbone and terminal declarations. It is portable
to native and WASM targets and has optional JSON schemas. This batch supplies the
scientific core, imported-annotation review, binary worker, native MCP access and
workspace controls with exact-event declaration storage. Graphical/platform
acceptance remains outstanding; this does not complete F34.

## Representation and scope

A `Declaration` binds an exact immutable `Molecule` revision and includes its
current declaration version, coverage, termini, position-indexed residue and
linkage records, and documentary provenance. No label or citation grants measured
chemical-state or cloning authority. Provenance remains user-supplied text.

`Coverage::Partial` leaves quantities unavailable. `ClosedModel` is an explicit
hypothesis that unlisted residues have canonical beta-D ribose chemistry and
unlisted bonds are 3′–5′ phosphodiesters. There is no default closed assumption.
Neither this declaration nor a formula proves the actual chemistry of a sample.

| Layer | Supported declarations |
| --- | --- |
| Nucleoside | Canonical, N6-methyladenosine, 5-methylcytidine, 5-methyluridine, pseudouridine, inosine, explicitly unknown identity |
| Sugar | Ribose, 2′-O-methyl, 2′-deoxy, 2′-fluoro, explicitly unknown chemistry |
| Backbone | Phosphodiester, phosphorothioate with Rp/Sp/unspecified stereochemistry, explicitly unknown linkage |
| Linear 5′ end | Hydroxyl, monophosphate, triphosphate, explicitly unknown chemistry |
| Linear 3′ end | Hydroxyl, monophosphate, 2′,3′ cyclic phosphate, explicitly unknown chemistry |
| Circular | An explicit closing 3′–5′ bond, with no linear terminal groups |

Residue indices are zero-based in the stored reference. A linkage's `after` index
identifies the upstream residue; the last position is valid only for a circular
seam. Known nucleoside profiles require their exact canonical parent: A for m6A
or inosine, C for m5C, U for m5U or pseudouridine. Inosine's A projection does not
infer an editing mechanism or its effects on translation or pairing.

Base and sugar declarations share one residue record, allowing an explicit
combination without overlapping duplicate entries. Methylation of the base and
2′-O-methylation remain distinct even when elemental formulas coincide.
Unknown names are never guessed into known profiles. Caps, other modifications,
cross-links, nonstandard connectivity, ion/adduct state, isotope enrichment and
chemical combinations outside the supported model remain unknown or unsupported.

## Quantities and identity

For a complete supported declaration, the analyzer computes integer CHNOPS/F atom
counts, a Hill formula, monoisotopic mass and a conventional average mass. It models
the fully protonated neutral-acid species, without counterions, hydration or mass
spectrometric adducts. It does not predict physiological charge, Tm, folding,
stability, enzymatic behavior or translation through modified residues.

The algorithm sums nucleosides and adds phosphoric acid minus two waters per
3′–5′ linkage. A terminal monophosphate instead adds phosphoric acid minus one
water. The triphosphate increment accounts for its ester and two anhydride bonds;
a terminal 2′,3′ cyclic phosphate consumes both terminal hydroxyls. A known blocked
2′ position therefore contradicts a declared 2′,3′ cyclic phosphate and refuses.
Unknown 2′ chemistry withholds quantities rather than proving that cyclization is
possible. Circular chains add their explicit seam and have no linear-end additions.

Phosphorothioate swaps one non-bridging oxygen for sulfur. Rp and Sp have the same
formula but different declared identities. Unspecified stereochemistry can leave
the formula complete while chemical identity remains unresolved; no mixture
probability or experimental product distribution is invented.

Partial coverage, ambiguous reference bases, unknown residue/sugar/linkage/ends
and retained unreviewed source modification evidence make aggregate formula and
mass unavailable. Reasons include counts and at most 100 positional examples;
full declarations remain intact. Known modification annotations and metadata are
listed as retained evidence, not silently erased by an empty closed declaration.
Review/mapping of those imported records is required for quantitative use.

## Imported annotation review

`annotation_candidates` inspects retained modification features without changing
an imported molecule. Supported exact, single-residue stored-strand
`modified_base` features use the [INSDC feature table, version 11.4](https://www.insdc.org/submitting-standards/feature-table/)
`mod_base` codes: `m6a`, `m5c`, `m5u`, `p`, `i`, `cm`, `gm`, `um`, `fm` and `tm`.
Normalization of external code case/whitespace preserves the original qualifier.
Unsupported codes (including `OTHER`), other feature kinds, uncertain/compound or
complement locations, features flagged after editing and mismatched reference
parents remain unresolved. No unsupported label is guessed into a known profile.

Each explicit `FeatureReview` references the feature index, its exact content
fingerprint and a residue entry in the declaration. A review succeeds only when
the declared position and every chemistry component specified by that annotation
agree. For example, `fm` requires both pseudouridine and 2′-O-methyl sugar; `cm`
constrains the sugar while permitting an explicitly declared additional base
modification. Source edits, changed feature contents, missing entries, conflicting
chemistry and repeated reviews refuse. Reviews never imply closed coverage.

Reports retain all source evidence and separately identify the unreviewed subset.
One reviewed feature cannot clear another feature or opaque modification metadata;
the latter remains unresolved. Review records are documentary constraints, not
experimental confirmation. At most 4,096 feature reviews are accepted, within the
declaration byte bound. The complete imported source remains independently bound.

## Binding and limits

The source hash matches the shared molecule revision, using a bounded streaming
serializer. Declaration, catalog and request hashes bind the complete source,
provenance, selected chemistry, canonical/modified formulas and mass constants.
Changing a sequence or source annotation invalidates its old declaration. Saving
declarations alongside exact history events is the intended integration boundary;
embedding a self-referential full source hash into its own molecule is avoided.

The current declaration version is 2 and analysis engine is `rna-chemistry-2`.
Unsupported versions refuse directly. `Report::validate_for` recomputes all results
from independently supplied source and declaration; a saved hash/report alone is
not a chemical-state or material certificate.

Bounds include one million reference bases, 100,000 residue records, 100,000
linkage records and 8 MiB of serialized declaration. The shared borrowed model
budget additionally bounds source text, annotations, qualifiers and location
depth; the streamed source serialization is capped at 32 MiB. Duplicate residue
or linkage declarations, incorrect parents, invalid coordinates, terminal linkage
indices and topology/terminus contradictions fail before any complete-model claim.

## Worker and MCP

Current worker version 116 adds `AnalyzeRnaChemistry` (job 193) and
`RnaChemistry` (result 183). Both stateless and retained-session dispatch validate
the complete source/declaration context before analysis. Requests and responses
use the ordinary bounded binary envelope, without widening workspace admission.
Response decoding checks the current catalog/context, declaration hash, bounded
shape, review-evidence subset and formula/mass consistency. These checks cannot
reconstruct source bases absent from a response: saved-record acceptance still
requires fresh `validate_for` against the independently retained source.

Native MCP exposes:

| Tool | Input and behavior |
| --- | --- |
| `list_rna_chemistry_profiles` | Current versions, catalog, supported nucleosides/annotation mappings, mass model and bounds |
| `inspect_rna_annotations` | Exact immutable `document` reference; candidate feature fingerprints and residue suggestions, with full JSON evidence |
| `analyze_rna_chemistry` | Exact immutable `document` plus typed current `declaration`; fresh native report, quantities or explicit unresolved issues |

The MCP schema requires coverage and termini rather than supplying a closed model.
Source bounds are checked before cloning a session document. Annotation previews
contain at most 100 rows and 80 characters of each raw code; full artifacts retain
the original strings. Analysis previews bound stereocenter examples and issue
positions. Publication obeys the existing session artifact/storage budgets. No tool
mutates a source, saves a declaration or installs a candidate/material. Thin WebMCP
and filesystem permissions are unchanged.

## Workspace controls

Open an explicit RNA document and choose **Biology → RNA chemistry**. The panel
pins an exact retained document event. **Use active RNA version** starts a new
partial declaration for the selected event; it never silently retargets a saved
hypothesis to equal sequence content elsewhere. Linear termini initially remain
unknown. Enter a model statement, choose coverage and termini, and add any residue
or linkage declarations. Positions shown in the form are one-based; native storage
remains zero-based. A bond after the final residue exists only on a circle.

Residue controls combine base and sugar chemistry with documentary provenance;
backbone controls distinguish Rp, Sp and unresolved phosphorothioate stereochemistry.
Unknown descriptions remain unknown, and a base/sugar combination is checked by
the native analyzer. Existing entries can be edited or removed. Imported annotation
reviews merge only their positive component constraints, allowing e.g. separate
m5C and 2′-O-methyl annotations to describe the same residue. Removing a residue
also removes its associated reviews and updates surviving entry references.

**Calculate RNA composition** runs through the existing isolated worker. Changes
to the pinned source, declaration or workspace invalidate the previous response;
cancellation closes the worker. The report shows formula/mass availability,
identity limits, unresolved source records and bounded positional diagnostics.
Report download retains the full declaration/evidence. Residue, linkage and
annotation lists show at most 50 rows per page.

## Exact-event declaration storage and sharing

**Save RNA declaration** records a named declaration after a successful native
analysis, including analyses whose quantities remain unavailable. The current
`rna-declaration-v1` attachment contains record version 1, name, exact source event
and the current version-2 declaration. No report, mass, freshness flag or measured
state is cached. Opening restores the hypothesis against the independently retained
historical source; explicitly calculate to obtain quantities again. Saving and
deleting change only attachments, leaving source sequences and history intact.

Workspace archive validation rejects dangling/foreign event bindings, changed
source revisions, unsupported record/declaration versions, incompatible annotation
reviews, extra fields, invalid roles/hashes and oversized records. Serialization is
bounded to the declaration limit plus 4 KiB for its event/name wrapper. The ordinary
workspace attachment and expanded/archive limits still apply. Original annotation
metadata is preserved; no sequence/material certificate is created.

Whole-workspace saves, browser autosave and offline desktop storage retain these
attachments. A single-document `.cg2`, GenBank or FASTA export contains the reference
document, not these separately named hypotheses; use a workspace file or native
selected-library package to retain them. Selected-document sharing uses closure engine
`library-document-closure-4`: declarations follow only their exact retained source
events, including retained ancestry. A declaration for another document with equal
content is omitted. There is no automatic application of a declaration to a newly
edited sequence and no compatibility reader for previous internal versions.

## Primary sources and fixtures

Canonical nucleoside formulas follow ChEBI's neutral
[adenosine](https://www.ebi.ac.uk/chebi/CHEBI:16335),
[cytidine](https://www.ebi.ac.uk/chebi/CHEBI:17562),
[guanosine](https://www.ebi.ac.uk/chebi/CHEBI:16750) and
[uridine](https://www.ebi.ac.uk/chebi/CHEBI:16704) records. Modified profiles retain
references to [m6A](https://www.ebi.ac.uk/chebi/CHEBI:21891),
[m5C](https://www.ebi.ac.uk/chebi/CHEBI:20607),
[m5U](https://www.ebi.ac.uk/chebi/CHEBI:45996),
[pseudouridine](https://www.ebi.ac.uk/chebi/CHEBI:17802) and
[inosine](https://www.ebi.ac.uk/chebi/CHEBI:17596).

Sugar substitutions are checked against
[2′-O-methylcytidine](https://www.ebi.ac.uk/chebi/CHEBI:19228) and
[2′-deoxy-2′-fluorocytidine](https://pubchem.ncbi.nlm.nih.gov/compound/101507).
The oxygen/sulfur substitution is documented by
[IDT](https://ensemble.idtdna.com/site/Catalog/Modifications/Category/8).
Mass constants use the listed common-isotope masses from
[NIST](https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl?ver=html)
and explicitly declared conventional average weights; natural-isotope variation
is not converted into a claim of exact average mass.
Terminal fixtures also use the neutral [AMP](https://www.ebi.ac.uk/chebi/CHEBI:16027)
and [ATP](https://www.ebi.ac.uk/chebi/CHEBI:15422) records.

Authored fixtures compare published monomer quantities, AMP/ATP formulas and
an independent nucleotide-condensation oracle across all 340 canonical sequences
of lengths 1–4 in both linear and circular form. Other fixtures check isobaric
identities, sugar/backbone substitutions, terminal/seam chemistry, unknown and
retained evidence, budgets, report/source tampering and exact-feature annotation review. Native/WASM checks do not
establish UI, vendor-notation or experimental acceptance.
