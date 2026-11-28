# Geneious Biologics antibody workflow and GUI evidence

Research checkpoint 2026-11-28. **Geneious Biologics supports substantial antibody
analysis**; it is a relevant workflow reference for Collagene F34. This audit uses
Biologics official support articles and eight downloaded official screenshots,
actually inspected. It does not attribute Biologics capabilities to Geneious Prime.
Screenshots document particular releases, not a live trial of the latest product.
No Collagene browser, application HTTP, screenshot, build or test checks were run.

## Visually observed controls and sequence presentation

| Official screenshot inspected | Direct visual observations |
| --- | --- |
| [Annotation menu](https://help.geneiousbiologics.com/hc/article_attachments/360059680612) | Folder tree on the left; breadcrumb and file list to its right. Top-level Export, Pre-processing, Annotation and Post-processing menus. Antibody annotation is a menu operation on selected files, not a permanent row of specialized sequence-view buttons. |
| [Numbering options](https://help.geneiousbiologics.com/hc/article_attachments/10388186364308) | Expanded Analysis Options groups numbering, variant annotation, protein statistics and liabilities. A numbering dropdown and separate number-display checkbox are adjacent. The visible AHo item explicitly uses **Kabat CDR definitions**; numbering and CDR boundaries cannot be treated as synonymous. |
| [Analysis options](https://help.geneiousbiologics.com/hc/article_attachments/20158105577876) | One compact scheme selector; checkboxes reveal optional functions. Liabilities occupy an editable configuration area; additional-feature settings are below. These are analysis-dialog controls. |
| [Numbered sequence](https://help.geneiousbiologics.com/hc/article_attachments/11785653015700) | DNA and translated amino acids share coordinates. Separate VDJ, FR, numbering, gene and variant tracks align vertically. Number labels are a pale blue track; FR is green. Right-side accordions group translation, annotation-type visibility, graphs, sidebar and preferences. Hover exposes an AA variant's source interval, query/database codons and residues, numbering position and HGVS. |
| [Gene hover/details](https://help.geneiousbiologics.com/hc/article_attachments/11765567503892) | Table selection drives the sequence view below. A gene tooltip gives species/chain, source accession/allele, identity, coverage and mismatch details. CDRs are red, frameworks green; germline differences occupy separate purple tracks. Several alternative D/germline matches remain visible. |
| [Missing-region alignment](https://help.geneiousbiologics.com/hc/article_attachments/360061861811) | Names/tree at left, residues and consensus centrally, right accordion controls. The missing part has an explicit **region not found** annotation; it is not filled with invented sequence. Sequence Viewer and Sequence Logo are adjacent tabs; search and zoom are above the canvas. |
| [Alignment sidebar](https://help.geneiousbiologics.com/hc/article_attachments/19896618211220) | Names and selectable metadata/heatmap columns remain aligned with sequence rows. A right-side filter searches available metadata. Sequence Viewer, Sequence Logo, Tree View and Info are tabs; data columns do not become floating canvas buttons. Assay-data support is observed here, not proposed as new ELN scope. |
| [Alignment export menu](https://help.geneiousbiologics.com/hc/article_attachments/26394496139924) | A single Export dropdown offers document, image, table and selected-sequence extraction. Selected rows use checkboxes/highlighting. |

Local research copies are under `/tmp/collagene2-geneious-antibody-research/` with
matching names: annotation-menu, numbering-options, cdr-options, numbered-sequence,
gene-tooltip, region-alignments, alignment-metadata, export-menu. These are external
vendor documentation assets; they were not copied into the MIT product.

## Documented capability inventory

**Numbering and boundary definitions.** Biologics documents IMGT, Kabat, Chothia,
Martin, AHo and numbering inherited from a reference database. IMGT is the default.
It identifies Kabat 1991, Chothia/Al-Lazikani 1997 and Martin 2008 interpretations;
these names require precise versioned rules. Database-derived numbering has special
reference-format requirements; arbitrary custom references do not automatically
provide IMGT numbering.
[Official numbering schemes](https://help.geneiousbiologics.com/hc/en-us/articles/360044625812-Numbering-Schemes)

Small custom boundary changes use signed per-framework start/end offsets, with
heavy and light settings separated; larger deviations can use custom annotated
references. The screenshot's AHo/Kabat combination is an observed concrete reason
to keep our scheme and CDR-definition choices separate.
[Boundary adjustments](https://help.geneiousbiologics.com/hc/en-us/articles/360052017751-Adjusting-CDR-definitions)

**Species and references.** Provided germline databases include human/mouse Ig and
TCR plus alpaca Ig; human, mouse and alpaca Ig are described as fully supported.
Other species are offered on request/as-is or through custom references. Germline,
whole-variable-domain template, linker and additional-feature databases have
different roles. This does not establish universal species support.
[Reference database scope](https://help.geneiousbiologics.com/hc/en-us/articles/360045069531-Understanding-Reference-Databases)

**Detection and molecular arrangements.** Annotation identifies FR/CDRs and closest
V/D/J matches. Settings include automatic heavy-versus-light, explicit single-chain,
linked heavy/light scFv, associated separate chains, opposite-direction chains and
tandem heavy/VHH-like constructs. Association ambiguity can remain unpaired or
produce alternative heavy/kappa and heavy/lambda pairs. Equal gene matches can be
retained as a group, split by frequency or treated as unknown. Long-read mode can
disable numbering. Those are documented operation semantics, not proof that all
molecules receive a unique domain assignment.
[Antibody Annotator](https://help.geneiousbiologics.com/hc/en-us/articles/360045069571-Antibody-Annotator)

**Liabilities and sequence detail.** Optional screening includes stops/frameshifts,
annotation/sequence-quality problems, deamidation, isomerization, cleavage,
oxidation, glycosylation, cysteine patterns and selected hydrophobic motifs;
configurable assets and penalties are summarized in a score and per-region
annotations. This is vendor heuristic scoring, not a measurement we should copy
or relabel as probability.
[Liabilities](https://help.geneiousbiologics.com/hc/en-us/articles/5324606679060-Antibody-Sequence-Liabilities)

Numbered-position rules distinguish an unexpected residue from a missing numbered
position. Their default conservation rules are based on human Ig and are less
applicable to other species. Preserve this distinction if implementing our own
position-based checks with published, attributed rules.
[Positional liabilities](https://help.geneiousbiologics.com/hc/en-us/articles/11185576529940-Positional-Liabilities-based-on-Antibody-Numbering)

**Ambiguous/no-domain UX.** The result table separately reports chain identity,
annotation completeness, reading frame and stop status; completeness alone does
not imply an in-frame, stop-free product. Closest-gene ties can list both genes.
Per-region nucleotide/AA sequences, lengths, germline mismatches and HGVS have
separate columns. Positional failures and missing regions remain explicit.
[Result columns](https://help.geneiousbiologics.com/hc/en-us/articles/13875423818004-Antibody-Annotator-All-Sequences-Table-Columns)
No screenshot of a completely non-antibody/zero-domain result was located; its
exact empty-state layout is **not observed**. Do not infer it from the partial-region
example, or convert annotation failure into proof that a domain does not exist.

**Comparison and export.** Post-processing alignment accepts selected regions or
whole sequences; region-aware alignment keeps corresponding regions together.
Missing regions are annotated; sequences with none of the selected regions can
be omitted. Documented exports include Geneious/GenBank/FASTA, PNG, residue tables
(CSV/XLSX) and selected aligned regions or available parent sequences. MUSCLE/MAFFT
and tree engines are vendor capabilities; our external-runtime prohibition and
explicit MSA-port deferral still apply.
[Alignment and export](https://help.geneiousbiologics.com/hc/en-us/articles/360045069471-Sequence-Alignment)

Separate comparison documents can compare shared regions/clusters across multiple
annotation results, with filters, a reference dataset and count normalization.
This is distinct from ordinary alignment and does not justify adding organizational
reporting, inventory or assay/ELN machinery to Collagene.
[Cross-result comparisons](https://help.geneiousbiologics.com/hc/en-us/articles/360045069851-Comparing-Results-across-Multiple-Experiments)

## Mapping to the current Collagene implementation

This section is a source audit, **not an observation of the running GUI**.

| Existing implementation | Concrete gap / proposed next behavior |
| --- | --- |
| Biology inspector has an expandable antibody section and advanced settings; protein and focused-alignment-row analysis run in workers. | Preserve menu/inspector grouping. Put reference/chain/region choices together; expose the routine analyze action there. Use a shared report layout rather than adding independent buttons around the canvas. |
| Three schemes: IMGT, Kabat1991, Chothia1997. Separately selected regions: IMGT, Kabat1991, canonical-chart Chothia, North2011. | Martin and AHo are real alternative-scheme gaps. Implement exact positions and independent CDR definitions together with unambiguous labels; never add menu options that silently reuse another scheme. |
| Versioned human IGH/IGK/IGL V/J corpus; conditional candidate tracks, retained alternatives and unresolved assignments. | Human is current scope. Mouse and alpaca/VHH references, broader validated domain/chain assignment and deliberate tandem/scFv domain/linker representation are real gaps. No D-gene assignment, nonhuman identity or chain pairing should be inferred from the current human V/J model. |
| Sequence overlay supplies CDR colors, sparse position labels, exact residue hover and unresolved underlines. Alignment overlay follows the selected row/candidate. | Add a compact domain/candidate summary strip and consistent FR/CDR/numbering track toggles. Offer full label tracks at useful zoom; retain sparse rendering when crowded. Geneious's separate coordinate/numbering/variant tracks are a useful layout reference. |
| Results expose candidate selection, CDR summaries, per-residue assignment table, diagnostics and full JSON download. | Region-level copy/export and a selected candidate comparison table would be useful. Export exact source coordinates, chosen scheme/definition, ambiguity and source/model identity with sequences. Full generic alignment support does not imply region-aware antibody alignment. |
| Liability screening has independent concrete/possible motif counts and supplied-CDR scope; no experimental-risk probability. | Add region-aware positional/constant-cysteine checks only where supported numbering/context exists. Keep ambiguous, missing and unsupported positions distinct; show clickable residue hits beside the same domain summary. Broader rule/preset coverage is useful; vendor scores are not a calibration shortcut. |
| Zero retained candidates currently receives a numeric status; candidate overlays and rows depend on available tracks. | Add an explicit empty-state reason: no supported candidate, unsupported species/model, partial evidence, ambiguous evidence or search truncation. Keep the source visible, and keep absence claims scoped to searched models/limits. A partial candidate must remain inspectable. |

Source locations: `crates/web/src/biology/antibody/{antibody_view,alignment_antibody_view,antibody_overlay}.rs`,
`crates/web/src/viewer/sequence/sequence_view.rs`,
`crates/antibody-numbering/src/{domain_tracks,v_mapping,references}.rs`,
[our numbering contract](../ANTIBODY_NUMBERING.md) and
[liability contract](../PROTEIN_LIABILITIES.md). The liability document has older
introductory wording about automatic numbering; the newer numbering implementation
and source show the actual conditional human-reference capabilities above.

Implementation priorities inferred from this comparison: improve coherent
candidate/domain/region presentation and empty states; add useful region extraction
and comparison; then broaden versioned scheme/reference/rule coverage. These are
product tasks, not a new broad calibration or acceptance project. Retain the user's
existing deferrals for external MSA tools and high-throughput BAM work.
