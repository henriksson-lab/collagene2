# Protein properties, scientific profile 1

`protein_properties::analyze_protein(&Molecule, &ProteinOptions)` is a pure Rust, native/WASM calculation for an explicitly typed, linear Protein document. No installed scientific program or online calculator is invoked. The report binds the complete molecule revision, exact stored sequence hash, model identifiers and all options. It reports source length separately from analyzed residue count. These are sequence/model estimates, not evidence of activity, expression, folding, solubility, mature processing or experimental identity.

## Chain and evidence assumptions

The default permits exactly one terminal `*`, explicitly records its removal, and rejects internal or repeated stops. The reject-stop policy also rejects a terminal stop. A stop-only document is not a peptide. No downstream sequence is discarded and no separate chains are joined automatically. Each modeled chain has free N/C termini; multimers, circular peptides, terminal blocking and unspecified modifications require other chemical models.

Composition retains exact symbol counts and fractions, including unresolved identities. `B` means D/N, `Z` E/Q and `J` I/L. `X` never receives an invented average identity. The default makes B/Z/X mass unavailable; selecting `CanonicalBounds` explicitly treats each X as one of the 20 canonical residues. Bounds enclose possible masses and absorption values, without probabilities or a claim that their independent extrema are jointly attainable. J is isobaric in this mass model and does not change the constant-pKa charge calculation, but its identity remains unresolved. U/O remain visible in composition; quantitative properties are unavailable because this profile has no audited model for them.

`ModificationState::AssumeUnmodified` declares the raw-chain model except for the separately specified intrachain disulfides. `KnownOrUnknownModifications` makes quantitative results unavailable. Recognized retained modification/cross-link annotations (`MOD_RES`, modified_residue, cross-link/cross_link, bond, LIPID, CARBOHYD, DISULFID) and explicit modification metadata likewise prevent an unmodified-chain numerical claim. Arbitrary free-text notes are not a chemical parser. Imported raw reference sequence does not establish its mature chemical state.

## Mass and disulfides

Neutral mass is the sum of elemental residue masses plus one H2O for free linear termini. The fixed average/monoisotopic element dataset and canonical residue formulas follow the primary [University of Washington Proteomics Resource table](https://proteomicsresource.washington.edu/resources/data-analysis/masses/), reviewed 2026-11-24. It is a declared numerical convention, not a claim to equal every calculator's atomic-weight revision. Average masses are natural-isotope convention values; monoisotopic masses use the dataset's light isotope convention. These are not charged-ion m/z values.

`Reduced` assumes no intrachain disulfides. `Count(n)` requires at least 2n explicit C residues; unknown identities cannot authorize a bond. `MaximumIntrachain` assumes floor(C/2) pairs, with an unpaired odd Cys; canonical X bounds widen the possible count conservatively. Every modeled bond removes 2H from mass, excludes its two Cys thiols from the charge model, and contributes 125 to the 280 nm extinction estimate. No interchain bond or actual residue pairing is inferred from a count. Retained explicit cross-link annotations are not silently reconciled with this ideal count model.

## Charge and theoretical pI

The implementation independently sums Henderson–Hasselbalch group occupancies. Basic groups contribute `n/(1+10^(pH-pKa))`; acidic groups contribute `-n/(1+10^(pKa-pH))`. The numerical dataset is the published default Epk values in the official [EMBOSS iep documentation](https://emboss.sourceforge.net/emboss/apps/iep.html): free N 8.6/C 3.6; Cys 8.5, Asp 3.9, Glu 4.1, His 6.5, Lys 10.8, Arg 12.5 and Tyr 10.1. Only these factual constants and the documented model were consulted; no implementation was copied or invoked. In this application J retains I/L semantics, not the older dataset's other residue-letter conventions.

A 64-step bracketed bisection finds zero model charge within pH 0–14. If there is no sign-changing bracket, pI is unavailable; a boundary is never substituted as a pI. Net charge at the selected pH can still be available. B/Z/X make charge and pI unavailable because their ionizable groups are unresolved. This is explicitly the Epk model, not a claim of ExPASy/Bjellqvist parity. Structural, ionic-strength, temperature and residue-context effects are omitted. A numerical pI is a theoretical model result.

## A280 estimate

The empirical composition formula uses `5500·W +1490·Y +125·cystine`, in L mol^-1 cm^-1 at 280 nm for the declared cysteine state. The primary [ExPASy documentation](https://web.expasy.org/protparam/protparam-doc.html) identifies the Pace calibration and its limitations. Absorbance at 1 mg/mL and 1 cm path is epsilon divided by average neutral chain mass. The model assumes no other chromophores; a zero coefficient is not proof of no actual UV absorption. Folding, conjugates and protein environment can change measured absorption, particularly for proteins without Trp. The computed number is not an observed spectrum.

## Bounds and validation

The default residue limit is 1,000,000; the hard limit is 10,000,000. Exceeding the configured bound rejects the operation instead of reporting a partial molecular mass. Preflight bounds source evidence to 32 MiB of text plus conservative per-node overhead, at most 1,000,000 entries and 64 nested location levels before source validation/hash serialization. Output is constant sized: at most 26 composition entries, fixed numerical fields and bounded reasons/warnings. Charge solving uses composition counts, not repeated full-sequence scans.

`ProteinReport::validate()` checks strict structural integrity, finite ordered intervals, source/model bindings, composition totals and completeness constraints. Serialized reports do not gain scientific authority: clients needing correctness recompute from the bound source/options. `complete` means all requested quantities and identities are resolved under the declared model; it never means empirical validation.

Independent tests transcribe published residue masses, check free amino acids and peptide-water stoichiometry, compare model pIs with independently computed 50-digit Decimal Newton roots, enumerate literal resolutions of ambiguous/disulfide cases, and test stops, source/model changes, unsupported chemistry, malformed reports and resource bounds. No other cloning application is a correctness oracle.

## Web panel and artifacts

The Biology inspector includes **Protein properties** for explicit Protein documents. It exposes stop handling, ambiguity policy, cysteine pairing/count, modification declaration, charge pH and the residue budget. Changing a source or option clears the previous report immediately and terminates its pending worker; an explicit Cancel control and unmount cleanup do the same. The interface displays intervals with units, separate unavailable reasons, retained composition and model assumptions, with a complete source/options-bound JSON download. A point mass for J remains visibly identity-unresolved. No property report creates a verified cloning plan or edits a molecule.
