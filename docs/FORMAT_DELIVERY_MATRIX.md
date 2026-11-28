# Format delivery matrix

Source audit, 26 November 2026. “Implemented” here means code paths exist, not
that all dialects, every foreign field or GUI acceptance are complete. Read,
regenerated write and exact-original download are different capabilities. This
matrix supersedes older plan paragraphs that describe an entire format family
as either finished or absent. Recent format adapters and compatibility UI/MCP hooks are now registered by the
parent; focused batch fixtures and native/WASM compiler checks pass. Browser acceptance
remains deferred. Each row identifies actual direction/profile gaps rather than
app-wide promises.

The scientific model remains independent of adapters. Imported geometry, history,
checksums and format metadata are documentary; only fresh native reconstruction
can establish material/replay authority. None of these migration tasks includes
ELN/inventory, external alignment engines, OpenGel or further BAM/index work.

| Family/profile | Core read | Regenerated core write | Web/native frontend | MCP surface at audit | Concrete remaining work |
| --- | --- | --- | --- | --- | --- |
| GenBank/DDBJ and ApE-style qualifiers | Yes, multi-record; `formats`/`gb-io` | Yes, annotated export + notices | File/paste/URL import and File export | `import_document` text; native `.cg2` export | RNA T/U encoding and citation-source replay are corrected; foreign application qualifiers are not complete application databases |
| FASTA/plain DNA/RNA/protein | Yes | Yes, sequence-only | File/paste/URL import; File FASTA export | Generic text import; native artifact export | No annotations/history in FASTA; explicit alphabet/declaration review is now available in File/Paste and MCP, and contradictory literal labels cannot override actual export type |
| Native single document `.cg2` | Versioned parser | Versioned writer | File save/open; native durable export | `import_sequence_file` current NativeDocument profile; `export_document`; desktop catalog preload | Single molecule is not whole workspace/history closure |
| Native workspace/history `.collagene` ZIP and design/candidate envelopes | Yes, bounded schemas and immutable graph validation | Yes | Working autosave, workspace open/save, selected sharing, native catalog | Dedicated candidate/design/run surfaces; not a generic filesystem importer | Preserve exact schema/version capability floors; history closure and whole-workspace sharing remain separate intentions |
| SnapGene `.dna` / selected `.prot` | Bounded DNA reader and NEW strict uncompressed linear protein packet21 reader | NEW selected DNA/protein reference writer with strict/explicit-loss policy | File DNA/protein import, supported DNA URL import, original download, NEW guarded File SnapGene export | `import_sequence_file` explicit SnapGene profile; `export_snapgene` | Compressed data, broader packets/history/version fidelity and vendor-produced corpus/application acceptance remain unverified |
| ABI/ABIF, SCF traces | Yes, channels/calls/qualities/peaks | Exact original; NEW selected SCF3 writer from retained real channels/calls/qualities | Trace file import/view/original export; NEW reviewed SCF export menu | `import_sequence_file` explicit ABIF/SCF profiles; `export_trace_scf` with precision/loss review | Selected regenerated SCF3 writer has native fixtures and MCP exchange; ABIF scores are not converted to four-channel SCF accuracies, and no signal is invented from sequence |
| Original Collagene XML | Yes, independent migration profile plus NEW explicit lower-only reference projection | Exact original only | Generic file/paste, original download, explicit lower-only migration review | Generic XML import; NEW `import_compat_sequence` explicit reference policy | Regenerated legacy XML writer and broader schema corpus; lower-only source is reviewed reference correspondence, never reconstructed material |
| DNA Strider/Serial Cloner XDNA v0 | Yes, explicit v0 binary profile | NEW reviewed v0 DNA/linear-protein reference writer | Generic file import, original legacy download; NEW reviewed File XDNA export | `import_sequence_file` explicit XDNA v0 profile; `export_xdna` with loss review | Binary writer fixtures and worker-backed controls have focused validation. Names/unsupported annotations require explicit loss review. Other versions need separate identification; filenames do not prove format |
| EMBL/ENA nucleotide text | Yes | Yes, explicit field profile | Generic file import; explicit File EMBL export | Generic text import; NEW export_annotated_sequence explicit writer profile | Broader header/feature fidelity and representative application corpora remain work; direct MCP annotated export is now implemented |
| GenPept + UniProt annotated protein text | Yes, independent adapters | Yes, explicit source-ID/date/mass/CRC policies | Protein-aware import and compact File export dialog/original source download | Generic text detection imports; NEW export_annotated_sequence explicit protein profiles | Broader FT vocabularies/representative sources; never alias protein GenPept to DNA GenBank |
| GCG/Wisconsin SINGLE sequence | NEW bounded `gcg_sequence`, explicit alphabet required | NEW sequence-only writer with source-ID/date policy | NEW guarded compatibility import/export dialogs and File action | NEW `import_compat_sequence` and `export_gcg` | Focused core/MCP checks pass; browser acceptance deferred; profile is not MSF or native Lasergene `.seq` |
| Aligned FASTA, CLUSTAL, Stockholm | Yes, explicit alphabet/row metadata | Yes, declared loss profiles | Retained alignment workspace/selection/export | Typed import/read/export tools | No aligner executed; unsupported annotation semantics remain evidence |
| PHYLIP strict/relaxed sequential/interleaved, MSF, NEXUS | Yes, bounded explicit dialects | Yes, options and loss notices | Alignment import/export selectors | Typed alignment tools/options | PHYLIP truncation is explicit; NEXUS/GCG variants outside profiles refuse |
| PIR/NBRF, GDE, SELEX alignments | Yes, new registered adapters | Yes within profiles | Alignment selectors wired in recent batch | Alignment enum/tools wired in recent batch | Unsupported circular PIR, GDE offset profiles/annotation or SELEX codes require explicit errors/evidence, not inferred biology |
| Oligo/primer CSV/TSV and ordering sheets | Yes, explicit column mapping/units | Yes | Primer/oligo import/edit/export workflows | Typed import/export tools and plan-oligo exports | Vendor-specific sheets are separate profiles; no inventory/stocks system |
| FASTQ | Yes, explicit quality encoding | Yes, retained/projected quality constraints | File quality import/trim/original and derived export | Typed FASTQ tools | Sequence edits cannot retain stale per-call qualities silently |
| GFF3, GTF, BED annotations | Yes, explicit reference/coordinates | Yes, bounded profiles | Annotation import/review/export | Typed annotation tools | Exact reference binding mandatory; no automatic transcript assembly from GTF rows |
| SBOL 3.1 structural JSON-LD | Yes, graph + explicit projection | Yes, supported graph or explicit losses | Explicit component/sequence review dialog, exact original/expanded graph/native export | Five graph/inspection/projection/export tools; native SDK validation passes | Broader RDF classes/properties/profiles remain unsupported; selected documentary subset sharing is implemented and reviewed |
| PDB/mmCIF coordinates | `collagene-structures`, bounded source-preserving coordinate/connection profiles | Exact retained source bytes; no regenerated coordinate writer | Tools → Structures, binary worker, WebGL view and workspace retention | `import_structure`, `inspect_structure`; explicit text up to 4 MiB, exact source/full-report artifacts | Graphical/hardware acceptance remains; no folding, symmetry expansion or inferred bonds. [Profile](STRUCTURE_IMPORT.md) |
| HTML/PDF protocol and SVG map exports | Native writers | Yes | Workflow-specific guarded export; SVG File menu | Protocol/plan exports; artifact surfaces | Human documents are not sequence interchange or verified material; PNG map writer/menu was not established by this audit |
| gzip/ZIP | Adapter/container-specific bounded support | Native ZIP/workflow-specific outputs | Native workspace ZIP paths | Dedicated format/session artifacts | No universal archive importer or arbitrary nested decompression is implied |
| SAM/BAM/VCF/BCF/reference datasets | Existing bounded modules/workflows | Existing selected/profile outputs | Existing sessions | Existing typed tools | Frozen/deferred area for this work; NO further BAM/BAI/CSI/CRAM/index tasks are proposed |

Core evidence: [generic sequence dispatch](../crates/core/src/formats.rs),
[legacy XML/XDNA](../crates/core/src/legacy_formats.rs),
[SnapGene](../crates/core/src/snapgene.rs), [traces](../crates/core/src/traces.rs),
[EMBL](../crates/core/src/embl.rs), [protein](../crates/core/src/protein_formats.rs),
[GCG](../crates/core/src/gcg_sequence.rs),
[SnapGene writer](../crates/core/src/snapgene_export.rs),
[alignment dispatch](../crates/core/src/alignment_formats.rs),
[extra alignment profiles](../crates/core/src/alignment_exchange_more.rs),
[annotations](../crates/core/src/annotation_formats.rs),
[oligos/FASTQ](../crates/core/src/oligo_formats.rs),
[SBOL](../crates/design-exchange/src/lib.rs).
Surface evidence: [App File/import flows](../crates/web/src/app.rs),
[alignment UI](../crates/web/src/alignment_view.rs),
[SBOL dialog](../crates/web/src/sbol_exchange.rs),
[MCP tool implementations](../crates/mcp/src/lib.rs).

## Practical cloning migration order

1. **Clone Manager**: prefer a reviewed Commented GenBank migration first. The
   [official CM12 manual](https://scied.com/downloads/CM12_Help.pdf) documents
   compact/commented GenBank, FASTA and batch conversion. It distinguishes native
   `.cm5` from these exports and describes additional notes/display information.
   Native CM5/CX5 readers are ABSENT here. Next independent deliverable: explicit
   Clone Manager GenBank profile/encoding fixtures (Windows text vs UTF-8, feature
   labels/colors/notes) over the existing parser. Then identify a pinned native
   version and independent sample corpus before designing a bounded binary reader.
2. **Vector NTI**: native MA4/PA4/OA4 archives/project databases are ABSENT here.
   Thermo Fisher provides an [official export route](https://www.thermofisher.com/us/en/home/life-science/cloning/vector-nti-software.html)
   and documents [GenBank/FASTA or protein exchange](https://www.thermofisher.com/order/catalog/product/A15339/faqs).
   Prioritize GenBank and protein-text migration fixtures with retained original
   bytes/unsupported qualifiers. Do not infer native archive support from being
   able to open one exported molecule; relationships and rich project data differ.
3. **DNASTAR/Lasergene**: native EditSeq/SeqBuilder `.seq`, `.mseq` and layout/project
   readers are ABSENT. DNASTAR's [official export documentation](https://www.dnastar.com/SeqMan_Pro_Help/Documents/commandsforsavingandexportingdata.htm)
   distinguishes native single/multi sequence, GenBank and GCG output. GCG is a
   useful independent bridge, but `.seq` is ambiguous: never select this adapter
   solely by that suffix. Next work: representative annotated GenBank and native
   version-signature fixtures, followed by a documented native profile if needed.
4. **Geneious, CLC, MacVector and other containers**: native database/project formats
   are ABSENT; prioritize their GenBank/FASTA/EMBL exports and original evidence.
   A safe intermediate export does not establish annotation/history completeness.
   Files with a shared extension must be identified by content/version/profile.
5. **Consolidate delivery surfaces**: provide one explicit compatibility import/
   export entry using typed format/alphabet/options and visible omissions, instead
   of adding one permanent toolbar button per format. Existing native/HTML/PDF
   history remains the fidelity escape hatch. Generic MCP sequence import still calls `formats::import_text`; the NEW
   `import_compat_sequence`, `export_snapgene` and `export_gcg` tools now provide
   explicit compatibility profiles and acknowledged output loss. Binary sequence/trace
   imports and XDNA/SCF exports now have direct MCP tools using the existing core
   adapters. Parser/UI/MCP and browser acceptance remain separate
   evidence, as documented in [compatibility controls](SEQUENCE_COMPATIBILITY_CONTROLS.md).

Native commercial-container work needs exact versions and independently sourced
fixtures. No third-party cloning application code was borrowed for this audit or
new GCG implementation. Further native proprietary writers should follow actual
migration demand rather than treating every parsed format as a writer obligation.

The selected XDNA/SCF writers pass fourteen focused native fixture groups
(eight XDNA, six SCF); annotated protein exchange passes nine and the MCP
annotated export tool passes two SDK groups. Worker-backed XDNA/SCF controls pass
their focused native UI-model and integrated web/worker WASM checks. This does
not establish browser acceptance, ABIF writing
or vendor application parity.
