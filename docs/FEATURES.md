# Features and supported formats

Collagene 2 is in active development. These tables summarize implemented capabilities; see the linked profiles for limitations and [implementation status](IMPLEMENTATION_STATUS.md) for remaining acceptance work.

## Workspace capabilities

- Sequence editing, annotations, translations and synchronized circular/linear maps.
- Restriction-site filtering, virtual digests, primer design and cloning workflows.
- Independent target verification, retained intermediate products and alternative cloning products.
- Saved history, primer orders and HTML/PDF protocols.
- MCP-assisted design, optional WebMCP, shared plasmid folders and offline desktop storage/sync.

See the [application guide](APPLICATION_GUIDE.md) for usage instructions, [desktop guide](TAURI_DESKTOP.md) for installation/build details and [server guide](SERVER_LIBRARY.md) for shared-library setup.

## Cloning methods

| Methods | Capabilities |
| --- | --- |
| Restriction/ligation | Sticky/blunt ends, ordered fragments and alternative products |
| PCR/mutagenesis | Added primer bases, mismatch/inverse PCR and sequential changes |
| Gibson/HiFi, In-Fusion | Overlap preparation and method-specific finishing assumptions |
| Overlap-extension PCR | Source PCRs, reciprocal extension, outer PCR and reusable products |
| Golden Gate | Type IIS preparation, fusion design and construct libraries |
| Gateway, TOPO | BP/LR, selected MultiSite profiles and explicit strand-transfer chemistry |
| TA/GC, end processing | Tails, phosphate treatment, blunting and oligo annealing |
| Domestication | Checked substitutions with inverse-PCR or overlap-extension preparation |

## File formats

| Formats | Support |
| --- | --- |
| GenBank/ApE, FASTA, EMBL | Import/export |
| GenPept, UniProt | Annotated protein import/export |
| SnapGene DNA/protein, XDNA v0, GCG | Selected import/export profiles |
| ABIF, SCF | Trace import, original download; SCF regeneration |
| Original Collagene XML | Import and original download |
| Aligned FASTA, CLUSTAL, Stockholm, PHYLIP, MSF, NEXUS, PIR, GDE, SELEX | Alignment import/export; no external aligner |
| FASTQ, GFF3/GTF/BED, primer CSV/TSV | Quality, annotation and oligo exchange |
| SBOL 3.1 JSON-LD | Selected structural profile |
| Native Collagene 2 | Documents, workspaces, designs and history |
| HTML/PDF, SVG | Protocol and map export |

[Format profiles and limitations](FORMAT_DELIVERY_MATRIX.md) · [NCBI/Addgene URL imports](URL_IMPORT.md)

## Platforms

| Platform | Delivery |
| --- | --- |
| Windows, macOS, Linux browsers | Leptos/WASM web application with WebGL2 |
| Windows, macOS, Linux desktop | Tauri; native CI builds, installed GUI acceptance in progress |
| Shared lab server | Rust service with authenticated plasmid folders |
| Mobile | Deferred |

