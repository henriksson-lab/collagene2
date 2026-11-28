# Translation and ORF model

`collagene_core::sequence::translation` is independently implemented MIT Rust. It performs no external calls and has no runtime dependency on an installed translator, aligner or design engine. The API is portable to native and browser WASM.

## Audited genetic-code profiles

The initial profiles are NCBI table 1 (standard), table 2 (vertebrate mitochondrial), and table 11 (bacterial, archaeal and plant plastid). Internal assignments and explicit accepted initiators come from [NCBI's genetic code tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi), whose displayed revision is September 2024. Tables 1 and 11 share internal assignments but have different initiation sets. Context-specific recoding, editing and programmed frameshifts are not inferred from these tables. Unknown table IDs return an explicit unsupported result.

DNA and RNA use the same codon model, treating U as the RNA equivalent of T. For IUPAC codons, every concrete expansion is checked. A unique shared amino acid is returned, otherwise `X`; a codon that could be either stop or sense is explicitly counted. Definite stops appear as `*`. These are conceptual translations, not expression predictions.

## Direct translation

`translate(&Sequence, &TranslationOptions)` returns `TranslationReport`. Options record the code, frame offset (0/1/2), strand orientation and initiation policy. Frame offsets count from the 5′ end of the selected orientation. Reverse orientation first reverse-complements the nucleotide sequence. Direct translation does not wrap at a circular origin or silently stop at the first stop; use explicit regions or the ORF API for bounded origin-crossing candidates.

- `Ordinary` reads every codon as an internal codon.
- `RecognizedStartAsMethionine` recodes the first codon as M only when all expansions are accepted initiators in that code.
- `RequireRecognizedStart` additionally rejects a missing, impossible or uncertain initiator.

Incomplete trailing nucleotides are counted and omitted, with a warning. Empty frame translations are valid reports unless initiation was required. Ambiguity and possible stops are separately counted, so a uniquely resolved synonymous codon is distinguishable from an unresolved protein residue. No padding or invented nucleotides are introduced.

## CDS translation

`translate_cds(&Molecule, &Feature)` uses the exact biological sequence specified by a CDS location, including compound and reverse/origin-spanning locations. It applies one phase offset after concatenating the biological parts, not a new frame at every exon. The extracted CDS is preflighted against the ten-million-nucleotide limit before allocation; repeated valid ranges cannot create an unbounded concatenation. Complement joins are traversed directly into one output buffer.

`codon_start` is interpreted as the one-based offset of the first complete codon, and `transl_table` selects the code, following the [INSDC feature table](https://www.insdc.org/submitting-standards/feature-table/). Missing values default to 1. Duplicate, flag-only, malformed or out-of-range scalar qualifiers fail explicitly. A phase greater than one is treated as an incomplete 5′ codon and does not trigger initiator recoding. A `partial` qualifier also suppresses complete-CDS claims.

The report contains the source revision and location, complete direct translation, CDS protein without a definite terminal stop, internal-stop count, warnings, and comparison with any stored `/translation`. Stored amino-acid annotations are never used as the computed answer. `complete_under_model` requires a definite initiator, complete codons, a terminal stop, no internal/possible stops or unresolved residues, no partial flag, and agreement with a stored translation if present. It is not experimental verification.

Fuzzy, ordered or unresolved locations, annotations already requiring review, and declared exceptions such as `transl_except`, `ribosomal_slippage`, frameshifts, pseudo/pseudogene or trans-splicing are unsupported. The ordinary sequence translator remains available for explicit exploratory reading frames without assigning those exceptional features an authoritative CDS interpretation.

## ORF candidates and coordinates

`find_orfs(&Molecule, &OrfOptions)` scans both orientations when requested. `AtgOnly` accepts definite ATG/AUG starts; `CodeStarts` accepts a codon only when every ambiguity expansion is an accepted initiator. Accepted alternative starts are represented as M in the ORF protein.

Each accepted start produces the candidate ending at its first definite stop or earlier uncertain-stop boundary. Nested starts are retained independently. ORF proteins omit the definite terminal `*`, while their nucleotide locations include its codon. Coordinates are zero-based half-open location trees in the original reference orientation. `frame` is the start's offset class in the selected strand, not an assertion of a global circular frame.

- `CompleteOnly` requires an accepted initiator and definite terminal stop.
- `IncludeLinearPartials` additionally returns linear-edge candidates, including a possible missing initiator at each of the three 5′ frame boundaries and an unclosed 3′ edge.
- `IncludeAllPartials` additionally returns accepted circular starts without a proven stop within one turn, or with an uncertain possible-stop boundary.

A possible-stop ambiguous codon forms an uncertain boundary; a candidate is not declared stop-free across it. A circular candidate consumes at most `floor(length / 3)` complete codons and never reuses a base through a second turn. This remains true when the molecule length is not divisible by three. A circular origin is not treated as evidence of a missing 5′ start. Partial flags and termination reason are explicit.

A rolling per-frame boundary scan is linear in source length and avoids a doubled genome allocation; it uses one bounded `u32` successor array. Codon outcomes are cached over IUPAC mask combinations. Results are limited to 10,000 ORFs and eight million materialized amino acids. Reaching a result budget sets `truncated` and warns that the result set is not exhaustive. Min-length and result-limit options are retained in the report. No minimum-length or sequence-content score is a biological expression guarantee.

## Validation and remaining scope

Targeted tests specify codons, phases, biological joins, reverse/origin coordinates, ambiguity, alternative starts and errors independently. An exhaustive collection of small circular canonical strings is also compared with a separate direct start-to-stop walk. No other cloning application is a correctness oracle.

Protein mass/pI, all remaining genetic codes, programmed recoding, exception-aware CDS repair, initiation-strength prediction and integrated fusion-frame verification remain future work. This module does not replace the deferred Primer3 or multiple-alignment engines.

## Web panel

The Biology panel sends explicit frame, selected-CDS and ORF requests to the Rust/WASM worker. Controls record initiation and partial-candidate policies rather than silently choosing them. Current molecule revisions, option generations and selected-CDS indexes gate returned results. Navigating away prevents an old asynchronous result from touching disposed component state.

Protein previews are bounded at 600 symbols and ORF lists at 100 visible entries; protein FASTA downloads retain the complete computed data. A candidate can be annotated at its full location as a predicted CDS requiring review; no stored protein translation is asserted by that action. Selecting an unannotated compound candidate displays its first reference span with an explicit notice, rather than pretending a single contiguous range is the whole ORF. Feature-directed translation uses the feature qualifiers, while manual frame controls apply to direct translation. Browser integration validation is separate from the core fixture suite.
