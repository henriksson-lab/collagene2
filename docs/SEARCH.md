# Unified document search

`search::search_document` is a readonly deterministic operation against an exact molecule revision. Search options choose `Sequence`, `Annotations`, `Primers`, `Enzymes`, or `All`. Results carry the full options, source revision, engine identity, bounded work/row counters and explicit completeness. Selecting a hit is a separate GUI action; searches never alter annotations or design verification.

Sequence queries use the document alphabet. DNA/RNA use IUPAC compatibility; ambiguity-compatible matches are distinguished from certain matches. RNA requires U and rejects T. Peptide search compares literal residues in the protein alphabet, forward only. Nucleotide both-strand search collapses palindromic queries to one result. Circular matches may cross the origin once, never repeat a molecule to satisfy an overlong query.

Annotation searches inspect feature kind and qualifier keys/values, including saved `primer_bind` annotations. Primer mode restricts those searches to saved primer annotations; it does not predict primer binding. Case-insensitive text matching uses Unicode lowercase. Whole-word matching requires non-alphanumeric/non-underscore boundaries. Sequence matching is always normalized by its alphabet and does not use text case/word settings. All mode skips an invalid sequence query with an explicit diagnostic and still searches annotation/enzyme names. Protein/RNA all-mode enzyme omission is explicit.

Enzyme name or recognition-string queries use the factual starter catalog, optionally restricted to explicitly named enzymes. Recognition matches retain their orientation/certainty and the actual top/bottom cut boundaries, including type IIS offsets. A cut outside a linear molecule is absent, not clamped; circular cuts wrap. Enzyme search requires DNA. This catalog is not a methylation, supplier-condition or experimental activity model.

Coordinates are zero-based half-open reference spans. Origin-spanning and compound feature parts retain biological order, including reversal of compound parts under `complement`. Fuzzy/order locations remain annotation evidence, never experimental confirmation. Unresolved features may match text without inventing coordinates. Feature indexes refer to the exact source revision. Stable hit IDs hash the source revision, engine identity, complete search options and typed hit identity.

Limits default to 1,000 hits, 100,000 feature rows and 10 million comparison/text-byte work units plus inspected feature, qualifier and location-node overhead; supported maxima are 10,000 hits, 100,000 rows and 50 million work units. Queries are at most 4,096 bytes. Reaching a work, row or hit bound returns explicit incomplete results across the requested sources. No silently complete first-source prefix is reported. Annotation locations with over 4,096 parts are omitted with incomplete status. The GUI paginates returned hits; MCP previews the first 100 and retains the full bounded report as a JSON artifact. Response decoding validates bounds, source binding and deterministic hit identities.

The current document search engine is `sequence-search-2`; superseded report
engine identities are rejected. Its work counter includes pattern copies/reverse
complements, singleton classification, KMP prefix initialization and symbol
comparisons. Concrete DNA/RNA uses one shared linear KMP kernel (T/U follow their
validated alphabet masks); protein sequence search remains literal. Overlapping
matches and every allowed one-turn circular start are retained. Compatible IUPAC
matching uses the same kernel's budgeted fallback, preserving the distinction
between possible overlap and definite subset matching. Exhausting work or output
capacity keeps the report explicitly incomplete.

The lower-level exact nucleotide API and all four exact-PCR binding searches
share that kernel. Exact search keeps its 100,000-hit cap and uses a fixed
200-million-work job limit; a PCR job shares that allowance across its four
orientation scans rather than resetting it for each primer. Those APIs cannot
return a partial vector as complete evidence: hitting either limit returns an
explicit error. Concrete late-mismatch inputs are linear, rather than being
rejected by a pessimistic `template_length * primer_length` preflight. Ambiguous
inputs that actually exceed the comparison budget must be refined. These are
embedded Rust algorithms, with no external aligner or primer program.
