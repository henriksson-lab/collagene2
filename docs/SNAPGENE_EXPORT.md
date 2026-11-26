# Selected SnapGene DNA/protein reference export

The independent MIT writer in `core::snapgene_export` serializes the **current**
native molecule. It is distinct from recovering an imported SnapGene file's exact
original bytes. It never copies opaque source chunks, stale sequence snapshots,
primer-binding caches or vendor cloning history into a changed document. Input
metadata/original source evidence remains unchanged.

## Public API

```rust
snapgene_export::review(&molecule) -> Result<Assessment, Error>
snapgene_export::export_snapgene(&molecule, Policy::Strict) -> Result<Report, Error>
snapgene_export::export_snapgene(&molecule, Policy::AllowOmissions) -> Result<Report, Error>
```

`Policy` is Copy/serde/JsonSchema. `Assessment` exposes `profile`, exact current
`source_revision`, `extension`, typed `losses` and `warnings`. `Report` adds `mime`
and binary `bytes`. All report/loss types are public and serializable. A loss has
`path`, `kind` and an explanatory `message`; kinds distinguish metadata/original
evidence, feature location/review/color/display and qualifier encoding. Strict
refuses any loss; inspect the typed assessment and explicitly choose AllowOmissions
before producing a lossy file. Mandatory safety/resource refusals apply to both
policies. Review does not serialize the sequence into a new output packet.

Profile `collagene-snapgene-reference-15-v1` writes a SnapGene cookie with format
revision15, uncompressed DNA type1/packet0 or protein type2/packet21, current name
in notes packet6 and selected annotations in packet10. All packet lengths are
big-endian u32. Linear/circular DNA has the observed double-stranded reference
property bit and topology bit; linear protein has a zero property byte. Returned
extensions are exactly `.dna` and `.prot`, MIME `application/octet-stream`. RNA and
circular protein are refused. A DNA duplex flag is a **reference representation**,
not a claim about experimental material, ends or modifications.

## Mapped annotation and note profile

Exact ranges and joins whose parts have one uniform orientation are exported with
one-based inclusive XML segment coordinates. Complement-of-join reverses biological
traversal; individually reversed parts are re-ordered explicitly before applying
a global reverse direction. Circular origin-crossing joins retain their actual
parts. Nested exact joins may flatten only while preserving source coordinates and
literal biological extraction order. Mixed-strand joins, order(), fuzzy endpoints,
unresolved/remote locations and reverse protein features cannot be flattened into
an asserted feature. Strict refuses; AllowOmissions omits the **whole feature** and
returns a typed loss.

Feature kind, label, duplicate qualifiers, flag qualifiers and literal values are
written as bounded XML. Attribute escapes include numeric tab/newline/carriage-return
references, so XML attribute whitespace normalization does not change values.
XML-invalid values and markup-like rich-text qualifier values are omitted only
with explicit losses; no HTML rendering or execution occurs in the writer.

A uniform valid `#RRGGBB` color is used without case conversion. Differing, missing
or invalid forward/reverse colors are never silently chosen/normalized: the segment
color is omitted with a loss and original color qualifiers remain literal.
Directionality0/3 hints are honored only when compatible with the actual current
location; conflicts/duplicate contradictions are losses, with literal hint
qualifiers preserved. The native annotation-review flag has no mapped field and
its loss is explicit. Display names may be derived from existing labels; importers
can add documentary display qualifiers, so exact native tree/qualifier-byte equality
is not claimed as a general round trip.

Native `definition`/`description`, `comment` and `accession` map to Description,
Comments and AccessionNumber notes. Competing values for one note are not merged
silently; the omitted field is a loss. Other metadata and original INSDC/SnapGene/
protein evidence are not embedded as current metadata. No old UUID, original source
confirmation claim or timestamps are guessed to be current identity. Experimental
confirmation is always explicitly0, never copied from producer/source claims.

## Chemistry and provenance boundary

Known native physical/duplex/modification, single/mixed-strand, TOPO, restriction,
variant and reference-only markers cause a hard refusal, including known generated
GenBank/EMBL producer markers. AllowOmissions cannot erase them and reimport an
apparently fresh intact duplex. Use the native document/workspace archive instead.
Known modified-protein metadata also refuses. Exact supported protein modification
annotations may remain documentary features; if exporting one would lose required
annotation evidence, both policies refuse rather than project an unmodified free
chain. Export does not establish native material, passing verification, kit
activity, methylation, host repair or experimental success.

The selected importer now also accepts the uncompressed type2/packet21 protein
profile. Reader and writer fixtures inspect independent literal packet bytes
before round-trip checks. Vendor-produced protein corpus and application
interoperability remain unverified; this profile does not imply broader packet
or chemical-state support.

## Bounds and confidence

Input/output sequence/file budget is32MiB; each XML packet is at most8MiB/100,000
nodes. Borrowed feature/model guards bound4096 features,16,384 qualifiers/location
nodes and32 nesting levels before recursive processing. Native metadata text is
bounded separately so retained original blobs are reviewed as omitted evidence,
not blindly decoded or forwarded. Source-revision hashing streams bounded JSON.
At most256 typed losses are accepted; excessive unsupported data refuses instead
of silently truncating the loss report. No partial output is returned on failure.

Layout facts were researched from the format investigators' own published notes:
[Damien Goutte-Gattat's Binary Sequence Formats](https://incenp.org/dvlpt/docs/binary-sequence-formats/index.html)
and the [sgffp author's layout specification](https://merv1n34k.github.io/sgffp/format-spec).
The latter documents topology bit0, duplex bit1, protein sequence type2/packet21,
notes/feature XML and directionality/coordinate observations. Only layout prose
was consulted; no application/parser implementation was copied, translated or used
as a scientific oracle. [SnapGene's developer-facing conversion page](https://www.snapgene.com/features/convert-file-formats)
is not a complete public binary specification. This selected profile is therefore
explicitly experimental interoperability, not vendor certification or an assertion
that every SnapGene release reopens every field.

Seven independent literal fixture groups cover exact cookie/flags/notes, reverse
origin and multipart biological ordering, XML escapes/duplicate/flag qualifiers,
protein bytes and alphabet/topology refusal, explicit color/display/location/
metadata losses, chemistry/modified-protein refusal, unchanged source recovery
following edits with no old opaque chunks, and model/XML bounds. No vendor GUI,
browser, external executable, compiler or test runner is invoked by this source
worker. The parent registers the module and runs coordinated native checks.
