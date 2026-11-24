# Basic primer tools

`collagene_core::primers` independently implements DNA oligo analysis, reference binding checks and anchored endpoint pair selection. No Primer3 source, port, thermodynamic library or external executable is used. Advanced folding/design engines remain deferred.

## Calculations and conditions

The optional perfect-duplex Tm uses factual enthalpy/entropy parameters from [SantaLucia 1998, Table 2](https://pmc.ncbi.nlm.nih.gov/articles/PMC19045/), terminal-base initiation and the self-complementarity entropy correction. The sodium entropy adjustment is `0.368 × (length−1) × ln([Na+])`. Tm is `1000 ΔH / (corrected ΔS + 1.987 ln(C)) − 273.15`; `C` is total strand concentration for a self-complementary oligo or total/4 for equal non-self-complementary strands. Conditions explicitly record sodium and total strand concentrations in mol/litre. Defaults are 0.05 M sodium and 500 nM total strands.

This module applies that model only to unambiguous, perfectly complementary DNA. Ambiguous inputs receive GC bounds and an unavailable-Tm explanation. No magnesium, dNTP, mismatch, modification or additive correction is inferred. The reported Tm is not a recommended PCR annealing temperature. Caller conditions are validated; sodium outside 0.01–1 M is rejected.

Annealing-domain Tm excludes the 5′ tail. Ordered-oligo GC, homopolymers and complementarity screens include it. Hairpin, homodimer and heterodimer screens enumerate ungapped contiguous antiparallel compatible runs; hairpins require three unpaired loop bases. They return exact intervals and 3′-terminal involvement. Ambiguous-compatible runs remain possible. These screens do not predict folding energies, gapped structures, equilibria or extension probability.

## Binding checks

`find_binding_sites` scans both primer orientations, including a circular origin, and keeps tails separate. Hits include source coordinates, a full `Location`, primer-oriented mismatch/ambiguous positions and the number of consecutive certain 3′ matches. Self-complementary primers retain both physical orientations. `BindingOptions` set allowed disjoint mismatches, a required certain 3′ run and work/hit bounds.

A compatible ambiguous position never becomes a certain exact base. Mismatch counts describe reference compatibility, not amplification. The report binds the source revision and options, and explicitly marks truncation: absence of more sites is not established after either budget is exhausted. Searches concern only this document, not a host genome or other laboratory templates.

## Candidate selection

`design_primers` takes an anchored half-open target interval. A circular interval uses an unwrapped end and covers at most one source turn. It enumerates requested forward and reverse annealing lengths, filters unambiguous candidates by GC, homopolymer, 3′ complementarity and pair Tm difference, then ranks by summed target-Tm distance plus pair Tm difference. Endpoint sequences are independently extracted through their full biological locations. Forward/reverse annealing domains may not overlap; the output length includes normalized tails.

Returned `PrimerCandidate` locations exclude tails, and reverse candidates use complement/join semantics where necessary. Saving them as `primer_bind` annotations should retain the ordered sequence, separate 5′ tail and calculation conditions as qualifiers. Candidate selection does not establish binding uniqueness: inspect binding evidence and use the independent PCR enumerator to review alternative products before use.

Oligos are limited to 200 ordered bases. Candidate annealing lengths are 8–80 with at most 64 tested lengths per endpoint; retained pairs are bounded to 100. Pair screening rejects workloads exceeding twenty million comparisons. Binding work is configurable up to 200 million base comparisons and 10,000 hits; defaults are twenty million and 1,000. These are application limits, not biological recommendations.

## Validation

`cargo test -p collagene-core --test primers --locked` includes manually summed dinucleotide/terminal/symmetry parameters, sodium/concentration checks, GC ambiguity bounds, exhaustive independent short-interval complementarity oracles, strand-oriented mismatch/3′ fixtures, circular origin binding, visibly truncated enumerations, deterministic ranked selection and independent PCR reconstruction of selected tailed primers. No other cloning application's output is a correctness oracle.

## Web workflow

The **Primers** inspector runs single-oligo analysis/binding and pair selection through binary Rust/WASM worker jobs. Changing the source or any calculation input clears obsolete results; request generation, source revision and component lifetime are checked before displaying a response. Reports download as JSON with source revisions and conditions. Binding-hit previews are limited to 100; downloads retain every returned hit and its truncation flag.

**Save primer annotation** records a `primer_bind` location, ordered oligo, separate annealing domain and 5′ tail, calculation notes and source revision. Ambiguous/mismatched binding annotations require review. **Load saved primer** uses the shared `primer_reuse` resolver pinned to the actual document revision and zero-based feature index. It checks all stored ordered/annealing/tail fields for agreement rather than treating the old binding location as current evidence. Generic imported primer annotations without oligo fields can supply their exact strand-oriented feature sequence only when geometry is supported, contiguous and unambiguous (including a single circular-origin crossing). Fuzzy, unresolved, clipped or reviewed geometry cannot invent a missing oligo. Deliberately stored mismatch primers remain literal oligos with warnings, not silently corrected to the reference. IUPAC bases remain ambiguous, and the ordered length is limited to 200 bases. The bounded saved-primer list exposes rejected entries and diagnostics, and reports omitted entries beyond its 1,000-item limit. The user must rerun analysis against the current document. This collection is document-based; a cross-document primer library and order-sheet formats remain separate delivery work.

The GenBank adapter explicitly recognizes valid DNA-valued `oligo_sequence`, `annealing_sequence` and `5prime_tail` fields and a valid 64-hex `source_revision`. It folds presentation whitespace inserted by annotated text wrapping only for those typed machine fields. Original records remain evidence; this profile does not normalize arbitrary unknown qualifiers or erase meaningful multiline notes.

Ranked candidate pairs can populate the PCR controls through **Use pair for PCR**. Analyzed/reused oligos can also be kept in forward/reverse slots and sent through **Use kept pair for PCR**. Those roles express user intent and do not constrain the independent enumerator to ignore another amplification orientation. Saving one candidate annotation changes the document revision and invalidates the old report; rerun selection before saving another. The independent PCR enumerator remains responsible for reconstructing all products under its stated binding model.


## Reuse from MCP

`resolve_primer` accepts `{document: {id, revision}, feature_index}` and returns `resolved: {source_revision, feature_index, primer, sequence_origin, warnings}`. The feature index is zero based in that exact immutable document. `primer` contains actual name, annealing domain and 5′ tail, suitable for an explicit `simulate_pcr` or `analyze_primer` request. This tool creates no document and makes no binding, thermodynamic or amplification claim. Historical `source_revision` qualifiers are retained evidence; they do not authorize current identity. Conflicting stored fields, unsupported inference geometry, stale revisions and non-primer features reject explicitly.

The source-authored SDK fixture resolves a literal tailed forward primer, infers a reverse primer from a literal complement interval, retains an intentional mismatch warning, rejects malformed/stale features, and passes the two resolved oligos into independent PCR reconstruction with a hand-specified product hash. These are authored fixtures; execution results are recorded separately.


The direct PCR inspector and exact-PCR rows in overlap and TA/GC workflows offer
saved-primer reuse across open documents. Select a source explicitly, inspect its
bounded primer list, and copy the desired primer into the forward or reverse
fields. Exact source event/revision checks run again on copy; names and separate
annealing/tail domains remain editable. Imported oligo-table documents use their
retained `primer_bind` annotations. A [small example](../fixtures/primer-reuse-demo.md)
shows the introduced forward tail and independently expected PCR sequence.

Direct PCR product previews now derive from the exact template event. Publication
stages all products and validates the combined workspace before changing open
documents. Template switches (including equal-content sources), edits, workspace
replacement and changed primer names/domains/tails invalidate obsolete responses.
This is retained reference-preview provenance, not scientific verification or
physical material authority. Browser acceptance remains deferred.
