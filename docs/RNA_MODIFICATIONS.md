# Declared RNA chemistry

`collagene-rna-chemistry` provides an independent MIT Rust analysis core for
source-bound RNA base, sugar, backbone and terminal declarations. It is portable
to native and WASM targets and has optional JSON schemas. This batch supplies the
scientific core; worker, UI, annotation review, workspace capture and MCP access
remain to be integrated. It does not complete F34 or change existing RNA documents.

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
Review/mapping of those imported records is still required for quantitative use.

## Binding and limits

The source hash matches the shared molecule revision, using a bounded streaming
serializer. Declaration, catalog and request hashes bind the complete source,
provenance, selected chemistry, canonical/modified formulas and mass constants.
Changing a sequence or source annotation invalidates its old declaration. Saving
declarations alongside exact history events is the intended integration boundary;
embedding a self-referential full source hash into its own molecule is avoided.

The current declaration version is 1 and analysis engine is `rna-chemistry-1`.
Unsupported versions refuse directly. `Report::validate_for` recomputes all results
from independently supplied source and declaration; a saved hash/report alone is
not a chemical-state or material certificate.

Bounds include one million reference bases, 100,000 residue records, 100,000
linkage records and 8 MiB of serialized declaration. The shared borrowed model
budget additionally bounds source text, annotations, qualifiers and location
depth; the streamed source serialization is capped at 32 MiB. Duplicate residue
or linkage declarations, incorrect parents, invalid coordinates, terminal linkage
indices and topology/terminus contradictions fail before any complete-model claim.

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
retained evidence, budgets and report/source tampering. Native/WASM checks do not
establish UI, vendor-notation or experimental acceptance.
