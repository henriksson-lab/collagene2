# Wisconsin/GCG single-sequence profile

The independent `core::gcg_sequence` adapter reads/writes a bounded single
ungapped sequence. It is distinct from GCG/MSF multiple-alignment support.
Original Wisconsin Package manuals describe the
[typed dividing line](https://bip.weizmann.ac.il/education/materials/gcg/using_sequences.html)
and [GetSeq output](https://bip.weizmann.ac.il/education/materials/gcg/getseq.html).
The latter supplies the independent eleven-call `GATTCCGATTG` / checksum `4920`
fixture. The [Biopython checksum API reference](https://biopython.org/docs/latest/api/Bio.SeqUtils.CheckSum.html)
provides a second independent uppercase/mixed-case checksum oracle (`5688` for
`ACGTACGTACGT`). Only documented values were consulted; no parser/library source
code was copied or linked.

Public API:

- `is_gcg(bytes)` admits only modern single-sequence magic or a bounded typed
  legacy dividing line. It does not choose an alphabet or validate all calls.
- `import_gcg(bytes,&ImportOptions { alphabet, fallback_name })` returns the
  ordinary `formats::ImportReport`. The caller MUST explicitly select DNA, RNA
  or protein. `Type: N` permits DNA/RNA but never guesses which; `Type: P`
  requires protein. Native alphabet validation rejects incompatible letters.
- `export_gcg(molecule,&ExportOptions { record_id, date_label })` returns an
  ordinary `formats::ExportReport`. Options are optional overrides of retained
  source fields; a synthetic record needs explicit ID and documentary date.
  No core clock, accession, placeholder date or biological event time is invented.
- `original_bytes(molecule)` returns the hash-checked exact source, not a
  regenerated representation of the current sequence.
- `checksum(sequence)` returns the Wisconsin uppercase-byte weighted checksum:
  cycling positions 1–57, sum modulo 10000. Its input must already be ungapped
  ASCII calls; it never strips arbitrary punctuation/annotations.

Import/Export options are serde and optional `json-schema` DTOs with unknown
fields refused. Modern `!!NA_SEQUENCE 1.0` and `!!AA_SEQUENCE 1.0` magic must agree
with the dividing line. Legacy no-magic input still requires explicit Type.
Length and Check (or documented Checksum spelling) are mandatory and compared
against actual calls. Optional numbered rows must name the exact next one-based
position. Alphabetic sequence groups and protein `*` stop calls are retained;
`*` is NOT a PIR record terminator here. Multiple records, multiple-alignment
magic, gaps, inline annotations and trailing positions refuse rather than being
silently stripped. Output uses numbered 50-call rows in 10-call groups.

Header/documentary text is not interpreted as native annotations or source
chemistry. Exact UTF-8 source bytes are retained as bounded base64 metadata chunks
with a BLAKE3 hash. Native revisions preserve that original through edits. Export
reports that features, primers, native metadata/attachments/history and original
bytes are omitted; circular topology and RNA's Type-N ambiguity have explicit
notices. Restrictive physical/reference-only material exports refuse instead of
laundering their ancestry through a format unable to carry those restrictions.
Known exact producer comments, if present in a source preamble, keep the existing
imported-chemistry restriction; no imported comment creates verification authority.

Limits: 4 MiB source, 1 MiB sequence, 64 KiB preamble/line, 100,000 body rows,
256-byte token ID and 128-byte documentary date label. The date label is bounded
source/caller text, not a validated timestamp or asserted biological provenance.

Six independent fixture groups cover the published checksum oracle, 57-position
wrap, explicit RNA/protein and real stop calls, length/checksum/type/coordinate
failures, unsupported annotation/alignment syntax, source retention/tamper,
synthetic export requirements and source limits. The implementation worker ran
formatting/diff checks only; parent compilation/fixture execution is pending.
The parent has registered generic detection (which refuses automatic alphabet
choice), guarded explicit-alphabet import/export dialogs, and MCP compatibility
import/GCG export tools. Their combined compiler/fixture batch is pending; no
browser acceptance is claimed. See [controls](SEQUENCE_COMPATIBILITY_CONTROLS.md).

The writer now calls the same comprehensive restrictive-evidence guard as the
selected SnapGene exporter. It refuses actual declared-source termini, imported
restriction/catalog evidence prefixes, retained original GenBank producer markers
and protein modification state as well as tracked/TOPO/variant/recombination
markers. The seventh fixture group constructs a real `DeclareSourceTermini` plan
and verifies its product before checking refusal; it also covers documentary
marker variants. This is NOT the stricter source-termini declaration preflight:
ordinary RNA, protein, circular reference topology and ambiguity calls remain
supported under the existing GCG omission policies.

For protein sources GCG additionally reuses
`protein_properties::has_modification_evidence`: MOD_RES/modified_residue,
DISULFID, cross-link/cross_link/bond, LIPID and CARBOHYD annotations cannot be
silently dropped and reimported as an unmodified free chain. The eighth literal
fixture covers all these feature types and a normal DOMAIN loss-control. This
restriction is specific to GCG's sequence-only output; it does not narrow the
SnapGene writer's separately reviewed faithful chemical-feature representation.
