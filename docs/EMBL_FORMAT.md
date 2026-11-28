# EMBL/ENA nucleotide adapter

The independently implemented Rust adapter lives in `crates/core/src/formats/sequence/embl.rs`. `import_embl` returns the shared import report; `export_embl` writes nucleotide records. It uses no runtime filesystem, network or external executable.

The framing, ID/SQ parsing, counters, header retention, evidence policy and writer are application code. Feature lines are adapted by their fixed columns into the existing `gb-io` standalone parser/writer because EMBL and GenBank share the INSDC feature/location grammar. Generated adapter GenBank text is discarded and never represented as original evidence. No cloning application's code was used.

## Supported behavior

- Modern seven-field IDs, explicit linear/circular topology, positive sequence versions, DNA and named RNA molecule types. Common four-field legacy IDs are accepted; an omitted legacy topology defaults to linear with a warning. Legacy SV headers map to the modern ID version; a regenerated legacy record uses a local default version of 1 if no version was supplied.
- Terminated single/multiple records, SQ nucleotide lengths, cumulative sequence-line counters and supplied nucleotide-composition counts. Sequences use validated IUPAC letters. Records with unknown/placeholder molecule types, protein units or unsupported ID layouts are rejected rather than guessed.
- INSDC compound/complement/order/fuzzy locations, unknown locations retained by the common parser, ordered duplicate qualifiers, flag qualifiers and multiline qualifier values. The parser retains multiline value boundaries as newlines. Exact feature extraction still rejects fuzzy/order/unresolved locations in the scientific model.
- DE maps to the common `definition` metadata. Other original headers, including references and unknown two-letter fields, remain in original evidence and are carried into regenerated EMBL exports unchanged. Unknown fields generate an import warning. Header values are not application instructions or permission grants.
- Each molecule retains only its own original record under `embl.original_record`, plus ID field metadata. The full multi-record input is not copied into every molecule. Original EMBL is deliberately not stored in the `original_genbank` field. Metadata evidence contributes to the normal native-file/workspace size budget.
- An unchanged molecule exports its original record exactly. Regeneration uses current name, topology, sequence, definition and features, while retaining the other original header lines. New molecules receive a local modern ID and a generated nucleotide-composition header. Feature review flags are rendered as explanatory qualifier notes by the shared writer.
- RNA is stored internally with U. For an RNA-declared record, supplied T is explicitly normalized to U with a warning; the original text preserves the supplied spelling. Regenerated nucleotide EMBL sequence data encodes U as T and counts it in the T field. This is an explicit adapter normalization, not alphabet guessing.

## Conservative limits and fidelity boundaries

Input/output are bounded at 32 MiB, individual lines at 1 MiB, records at 1,000, sequences at 10 million symbols, and feature/location structures by the existing core limits (including the 32-level location preflight). Record terminators and declared lengths must be consistent. Malformed input is rejected atomically.

CO assembly-only/reference-resolving records, unspecified ENA submission templates, protein/Swiss-Prot records and extended nonstandard ID/SQ layouts are unsupported. There is no remote accession resolution or ENA submission validation. The generated file is an exchange representation, not a claim that an edited accession constitutes a new database submission.

References, comments, assembly maps and unknown headers may contain coordinates. Regenerated export rejects a changed sequence when such headers are retained, until their remapping is implemented. It does not silently export those fields as if they had been remapped. For simple records whose retained headers are safe, an edited-sequence export includes a visible CC provenance note. Original evidence remains available through the native format.

The adapter retains headers in EMBL; it does not interpret every header as an editable database model. Export to a different format is a separate fidelity boundary: the common GenBank exporter does not map arbitrary EMBL headers, and generic native metadata/plan/chemistry fields have no blanket EMBL representation. Use `.cg2` to retain all application metadata and original evidence. An annotation/sequence round trip is not complete application-state round-trip support.

## Validation and sources

`cargo test -p collagene-core --test formats sequences::embl:: --locked` covers independently authored synthetic records, origin-spanning reverse features, fuzzy/order locations, duplicate/multiline/flag qualifiers, exact original export, regenerated header retention, malformed counts/terminators, RNA normalization, legacy IDs, edit/reference guards and excessive nesting. No other cloning application acts as a correctness oracle.

Primary format references:

- [ENA flat-file example](https://ena-docs.readthedocs.io/en/latest/submit/fileprep/flat-file-example.html) — field framing, feature and sequence blocks, reference positions and line counters.
- [ENA sequence-flatfile guidance](https://ena-docs.readthedocs.io/en/latest/submit/fileprep/sequence-flatfile.html) — annotated nucleotide records and INSDC feature conventions; submission templates contain placeholders that this archive-record adapter does not infer.
- [EMBL Nucleotide Sequence Database in 2006](https://academic.oup.com/nar/article-pdf/35/suppl_1/D16/3869521/gkl913.pdf) — modern ID field layout introduced with release 87.
- [INSDC feature table](https://www.insdc.org/submitting-standards/feature-table/) — common location and qualifier semantics.
