# Binary worker protocol

The browser and native worker share `core::protocol` DTOs and deterministic scientific engines. Request magic is `CGQ\0`; response magic is `CGR\0`. The current writer uses schema **67** and reads supported schemas 1–66. Candidate files write schema 15; workspace plan schema is 13; the verifier is science-18.

Each envelope has a 42-byte header: four-byte magic, little-endian `u16` version, little-endian `u32` payload length, 32-byte BLAKE3 payload checksum, then a Postcard payload. A request contains `id`, exact expected source `revision`, then a typed job. A response returns the same `id`/`revision` and a typed result or error. Compound/circular coordinates and scientific options are typed data, not JSON strings. Unknown versions, invalid lengths/checksums, trailing payload, stale source revisions and variants claiming a schema older than their introduction fail explicitly.

Ordinary payloads are limited to **32 MiB** during bounded serialization and before deserialization. Only native workspace save/load variants (schema 6+), the explicitly whitelisted cloning-run variants (schema 12+), digest installation (schema 20+), catalog digest installation (schema 21+), cloning output finalization (schema 39+), construct read capture/replay (schema 42+), read review (schema 43+) selected-batch begin/capture (schema57+), and atomic release/archive workspace cases (schema58+) can use the separate **128 MiB** message cap. A borrowed prefix probe checks those exact tags and bounded revision text before allowing a large payload; ordinary scientific jobs retain 32 MiB. The checksum provides corruption detection, not authorization.

## Stable job tags

Postcard enum indexes are append-only. Existing tags and supported historical layouts remain intact.

| Tag | Job | Minimum schema |
|---:|---|---:|
| 0 | Search (legacy nucleotide search) | 1 |
| 1 | Restriction | 1 |
| 2 | Pcr | 1 |
| 3 | Verify | 1 |
| 4 | DuplexDigest | 2 |
| 5 | Translate | 3 |
| 6 | TranslateCds | 3 |
| 7 | FindOrfs | 3 |
| 8 | Align | 3 |
| 9 | PcrDocuments | 4 |
| 10 | AnalyzePrimer | 4 |
| 11 | DesignPrimers | 4 |
| 12 | ExportProtocolPdf | 4 |
| 13 | GoldenGate | 5 |
| 14 | MismatchPcr | 5 |
| 15 | InversePcr | 5 |
| 16 | PredictGel | 5 |
| 17 | DiagnosticDigests | 5 |
| 18 | MismatchPcrPreview | 5 |
| 19 | SaveWorkspace | 6 |
| 20 | LoadWorkspace | 6 |
| 21 | ValidateTraces | 6 |
| 22 | ScanFeatureLibrary | 7 |
| 23 | ApplyFeatureLibrary | 7 |
| 24 | ImportAnnotations | 7 |
| 25 | ApplyAnnotations | 7 |
| 26 | ExportAnnotations | 7 |
| 27 | RepairPreview | 8 |
| 28 | AnnealPreview | 8 |
| 29 | AnnealedDocument | 8 |
| 30 | ViewTracks | 9 |
| 31 | SearchDocument | 9 |
| 32 | PrepareTracks | 9 |
| 33 | CachedViewTracks | 9 |
| 34 | ImportAlignment | 10 |
| 35 | ExportAlignment | 10 |
| 36 | ExtractAlignmentRows | 10 |
| 37 | SaveAlignment | 10 |
| 38 | LoadAlignment | 10 |
| 39 | ProteinProperties | 10 |
| 40 | TailPreview | 11 |
| 41 | CheckOverlapProfile | 11 |
| 42 | RecordCloningRun | 12 |
| 43 | ReplayCloningRun | 12 |
| 44 | DesignCodons | 13 |
| 45 | CheckCodonDesign | 13 |
| 46 | DotPlot | 14 |
| 47 | AssembleContigs | 15 |
| 48 | AssessMethylation | 16 |
| 49 | DesignCrisprGuides | 18 |
| 50 | SearchCrisprOffTargets | 18 |
| 51 | DesignEditTemplate | 18 |
| 52 | CheckEditTemplate | 18 |
| 53 | BeginGenomeFasta | 19 |
| 54 | PushGenomeFasta | 19 |
| 55 | FinishGenomeFasta | 19 |
| 56 | BeginGenomeNative | 19 |
| 57 | PushGenomeNativeBlock | 19 |
| 58 | FinishGenomeNative | 19 |
| 59 | AbortGenomeImport | 19 |
| 60 | RegisterGenome | 19 |
| 61 | GenomeIntervalGuides | 19 |
| 62 | StartGenomeSearch | 19 |
| 63 | AdvanceGenomeSearch | 19 |
| 64 | PageGenomeSearch | 19 |
| 65 | CancelGenomeJob | 19 |
| 66 | ReadGenomeInterval | 19 |
| 67 | InspectRestriction | 20 |
| 68 | MaterializeDigestFragment | 20 |
| 69 | InstallDigestFragment | 20 |
| 70 | InspectCatalogRestriction | 21 |
| 71 | MaterializeCatalogDigestFragment | 21 |
| 72 | InstallCatalogDigestFragment | 21 |
| 73 | InspectRecombination | 22 |
| 74 | CheckRecombinationTarget | 22 |
| 75 | SelectRecombinationComponent | 22 |
| 76 | ScanRestriction | 23 |
| 77 | PrepareReadAlignments | 24 |
| 78 | ReadAlignmentPage | 24 |
| 79 | CheckReadAlignment | 24 |
| 80 | ExportReadAlignments | 24 |
| 81 | DisposeReadAlignments | 24 |
| 82 | SimulateTopo | 26 |
| 83 | InspectTopoOrientations | 26 |
| 84 | ProjectTopoProduct | 26 |
| 85 | CheckTopoTarget | 26 |
| 86 | ProjectAnnotatedTopo | 27 |
| 87 | CheckAnnotatedTopoTarget | 27 |
| 88 | PrepareVariants | 29 |
| 89 | VariantPage | 29 |
| 90 | CheckVariantReference | 29 |
| 91 | PreviewVariantAllele | 29 |
| 92 | ExportVariants | 29 |
| 93 | DisposeVariants | 29 |
| 94 | PreviewOverlapRecipe | 31; catalog preparation requires 37 |
| 95 | PreviewLigationRecipe | 32; catalog preparation requires 37 |
| 96 | InspectAssemblyJunctions | 33; nested plan minimum also applies |
| 97 | DesignOverlapPrimers | 35 |
| 98 | BuildAssemblyIntent | 36 |
| 99 | SelectDigestFragment | 38 |
| 100 | ExportPlanOligos | 38 |
| 101 | ExportProtocolHtml | 38 |
| 102 | FinalizeCloningOutputs | 39 |
| 103 | DesignGoldenGatePrimers | 40 |
| 104 | BuildMutationIntent | 41 |
| 105 | CaptureReadCheck | 42 |
| 106 | ReplayReadCheck | 42 |
| 107 | PreviewReadReview | 43 |
| 108 | CommitReadReview | 43 |
| 109 | BuildEditedAssemblyIntent | 44 |
| 110 | DesignEditedAssemblyPrimers | 44 |
| 111 | DesignSelectiveAssemblyPrimers | 45 |
| 112 | DesignGoldenGateJunctionPrimers | 46 |
| 113 | InspectPreparedGoldenGatePool | 47 |
| 114 | SelectPreparedGoldenGateCandidate | 47 |
| 115 | InspectGoldenGateReaction | 49 |
| 116 | SelectGoldenGateReaction | 49 |
| 117 | PreparePcrOutputs | 50 |
| 118 | PrepareMismatchPcrOutputs | 51 |
| 119 | CompareSequences | 52 |
| 120 | CopyReferenceFragment | 53 |
| 121 | PasteReferenceFragment | 53 |
| 122 | TranslateReference | 54 |
| 123 | PrepareTranslatedReference | 54 |
| 124 | PredictPinnedGel | 55 |
| 125 | ExportPairwiseComparison | 56 |
| 126 | SearchPairwiseComparison | 56 |
| 127 | PreparePairwiseContext | 56 |
| 128 | BeginCombinatorialBatch | 57 |
| 129 | AdvanceCombinatorialBatch | 57 |
| 130 | CaptureCombinatorialBatch | 57 |
| 131 | BeginCombinatorialRelease | 58 |
| 132 | AdvanceCombinatorialRelease | 58 |
| 133 | InspectCombinatorialRelease | 58 |
| 134 | FinishCombinatorialRelease | 58 |
| 135 | ExportCombinatorialOutputs | 58 |

Result tags 0–29 retain their selectors (verification gains explicit versioned step evidence in schema 25) (`Search`, `Restriction`, `Error`, `Pcr`, `Verification`, `DuplexDigest`, `Translation`, `CdsTranslation`, `Orfs`, `Alignment`, `PcrDocuments`, `PrimerEvidence`, `PrimerDesign`, `ProtocolPdf`, `GoldenGate`, `MismatchPcr`, `InversePcr`, `Gel`, `DiagnosticDigest`, `MismatchPcrPreview`, `WorkspaceSaved`, `WorkspaceLoaded`, `TraceValidation`, `FeatureLibraryScan`, `AnnotatedDocument`, `AnnotationImport`, `AnnotationExport`, `RepairPreview`, `AnnealPreview`, `PhysicalPreviewDocument`). Schema 9 appends `ViewTracks` at tag 30, `DocumentSearch` at tag 31 and `TracksPrepared` at tag 32. New response variants cannot be decoded under older version labels.

Schemas 1/2 had verification reports without an engine identity; their actual historical response layouts migrate to the explicit `unrecorded/legacy-worker-v1-v2` marker. Old reports never acquire the current engine. Nested plans in Verify/PDF jobs are also checked against the schema that introduced their operations, preventing newer operations from hiding inside an older outer job.

## Schema 9 search and tracks

`SearchDocument { molecule, options }` returns `DocumentSearch(SearchReport)`. It binds source revision/length/feature count, full options and engine identity, and reports work/rows/completeness. Stable hit IDs cover the revision, engine identity, options and typed hit content. The decoder checks those identities, coordinate/index bounds, kind-specific fields and configured limits. See [search semantics](SEARCH.md).

`ViewTracks { molecule, window, options }` returns `ViewTracks(TrackReport)` for one-shot use. The report binds revision, engine/version, window and complete options, and retains work and residue counters, coordinate spans, cut positions and residue/codon mappings. Decoder validation checks report bounds and the response revision. See [track semantics](VIEW_TRACKS.md).

`AnalysisSession` supports persistent browser workers. `PrepareTracks { molecule }` validates an immutable source and builds its retained index, then returns `TracksPrepared(revision)`. Replacement is atomic: an invalid preparation leaves the prior source intact. `CachedViewTracks { window, options }` sends only viewport/options and the expected source revision. The session rejects missing, stale or replaced sources and returns a normal track report without cloning or transferring the whole document on each scroll. A stateless `process_request` returns an explicit session-required error for these two jobs. Cache lifetime belongs to the worker/session; it grants no document-editing authority.

The caller still checks request IDs, selected document identity and current revision before installing asynchronous results. Track/search results are readonly evidence; choosing a search hit or changing the selected region is a separate user action. Cancellation of a stale GUI request does not imply a new scientific verdict.

## Schema 10 exchange and protein properties

Alignment jobs exchange existing aligned FASTA, CLUSTAL and Stockholm files; no alignment engine is invoked. `ImportAlignment { bytes, format, alphabet, name }` and `LoadAlignment { bytes }` bind `request.revision` to BLAKE3 of the exact incoming bytes. `ExportAlignment { document, format }`, `ExtractAlignmentRows { document, row_ids }` and `SaveAlignment { document }` bind the exact canonical native JSON content identity from `protocol::alignment_revision`, including source evidence. Import/load return `AlignmentDocument`; export returns `AlignmentExport` with explicit losses; extraction returns `AlignmentRows`; save returns `AlignmentSaved`. Result tags 33–36 and job tags 34–38 require schema 10. Typed responses validate original-source integrity and document/molecule bounds. A saved alignment reply must additionally hash exactly to the canonical source revision carried by its response, preventing substitution of another valid alignment. These jobs retain the ordinary 32 MiB message cap, including envelope payload overhead; the inner attachment limit does not increase that transport budget.

`ProteinProperties { molecule, options }` at job tag 39 returns `ProteinProperties(ProteinReport)` at result tag 37. Both require schema 10. The source molecule revision must match the request, and the decoded report revision must match the response. Explicit chemistry/ambiguity assumptions, model identities, unavailable values and bounds remain typed. It does not update the sequence or confer verification authority. Those schema-10 jobs do not independently change candidate or verifier semantics; the later schema-11 method additions are described below.

## Schema 11 physical tails and overlap profiles

Jobs append `TailPreview` at tag 40 and `CheckOverlapProfile` at tag 41; results append `TailPreview` at tag 38 and `OverlapProfile` at tag 39. All require schema 11 and retain the ordinary 32 MiB message budget. Earlier tags and layouts are unchanged.

`TailPreview { molecule, recipe, base }` binds the source molecule revision. `TailSubstrate::Digest { enzymes, fragment_index, blunt_before }` reconstructs complete cleavage and optional explicit blunting. `ExactPcr { forward, reverse }` requires one exact product from the embedded competing-binding search. The result retains recipe, source revision, actual source/product duplexes, retained strand mappings and non-templated additions. Decoding replays the tail report and checks the preview's source/base/product bindings. A saved preview remains evidence; its metadata cannot authorize later physical transformations.

`CheckOverlapProfile { molecules, reverse, overlaps, circular, options }` binds the complete ordered list of original source revisions, orientation flags, homology lengths, topology and profile options using `overlap_profile_revision`. Reports retain original-reference strand coordinates. Exact endpoint matching is distinct from published advisory length guidance. Repeat-scan limits remain explicit, and even a complete scan does not establish exhaustive alternative assembly or experimental compatibility. The verifier recomputes checks from bound sources before sequence reconstruction.

Candidate files now write schema 5 for `Tail3Prime` and `ProfiledOverlap`; schemas 1–4 remain readable for operations they introduced. Native workspace plan schema 3 introduces those operations, retaining compatible plan schemas 1/2. The scientific verifier is science-7, including a blocked 3′ physical-substrate barrier for both legacy and profiled overlaps. New operations cannot be mislabeled as older candidates, Verify/PDF requests, inline workspace requests/results or native archive plans. Operation/orientation/recipe/primer/sequence deserialization rejects unsupported unknown fields rather than silently discarding requested chemistry.


## Schema 12 retained cloning runs

`RecordCloningRun { workspace, plan, target_binding, input_bindings, allow_incomplete }` appends job42 and returns `RunCommit` at result40. `ReplayCloningRun { workspace, run_id }` appends job43 and returns `RunReplay` at result41. All require schema12; tags/layouts1–11 remain unchanged. Candidate schema5, workspace plan schema3 and science-7 are unchanged because this feature records and replays existing operations.

Both request revisions are `run_lineage::workspace_fingerprint(&workspace)`, covering the entire immutable source workspace including historical revisions, event bindings, candidate data, attachments and navigation. The commit receipt and replay report contain `source_workspace_fingerprint`; decoding requires exact agreement with the response revision. Core validation checks retained run structure; recording always reconstructs afresh and replay never trusts a stored status as current authority.

Only request tags19/20 (minimum6) and42/43 (minimum12), and result tags20/21 (minimum6) and40/41 (minimum12), may exceed32MiB, bounded by128MiB. A borrowed id/revision/tag prefix checks the tag's minimum schema before allocating/deserializing any large body. Replay can carry bounded reconstructed products without lifting the ordinary `Verification` response cap. Browser one-shot worker wrappers and the worker input transfer apply the128MiB absolute ceiling, then the protocol's exact prefix gate; the persistent viewport cache keeps its ordinary 32 MiB ceiling.

## Codon design and independent candidate checking (schema 13)

Jobs append `DesignCodons { molecule, options }` at tag44 and `CheckCodonDesign { source, candidate, options }` at tag45. Results append `CodonDesign(DesignReport)` at tag42 and `CodonCheck(CheckReport)` at tag43. These variants require schema13 before their bodies are deserialized and retain the ordinary 32 MiB envelope limit; they are not workspace-bearing tags. Existing schemas1–12 retain their layouts and version floors. The request and response revision is the exact source molecular revision. Reports additionally bind the declared options revision and candidate revision. Callers reject stale option/source/candidate replies against their exact submitted snapshots. Native audit paths can use `DesignReport::validate_for_source` or `CheckReport::validate_for` to recompute the checker; the web uses `validate_bindings_for_source` / `validate_bindings_for` after its separate checking worker job, avoiding constraint scans on the main thread. Source evidence is bounded before revision hashing.

The codon design model has a separate engine identity. Candidate5, workspace plan3 and cloning verifier science-7 are unchanged: a designed synonymous sequence is not a new physical operation or a cloning verification result. Every returned design candidate passed the independent sequence/translation/constraint checker. A found candidate does not establish global optimality when search was pruned. [Scientific profile](CODON_DESIGN.md).


## Window dot plots (schema 14)

`DotPlot { molecule, query, options }` appends job46; `DotPlot(DotPlotReport)` appends result44. Their minimum schema is14, checked using the borrowed tag prefix before the variant body is deserialized. Both retain the ordinary 32 MiB payload limit. Existing schemas1–13 and all earlier tags/layouts retain their minimum-version rules. Candidate5, workspace plan3 and cloning science-7 are unchanged.

The envelope revision binds the exact reference molecule. `dotplot::source_revisions(reference, query, options)` bounds both sources and options before hashing; the report separately binds reference, query and options revisions. `validate_bindings_for` checks those identities against submitted immutable snapshots without rerunning the scientific scan. A stale query or changed options invalidates display/export even when the reference remains unchanged. Reports retain original-reference/query coordinates and compound circular spans, definite versus ambiguity-compatible match evidence, work/index/hit limits and stride sampling. A sampled or truncated plot cannot claim exhaustive full-grid coverage. It is window evidence, not an alignment path, homology proof or physical ancestry. [Scientific profile](DOTPLOT.md).


## Reference-free read contigs (schema 15)

`AssembleContigs { reads, options }` appends job 47; `Contigs(ContigReport)` appends result 45. Both require schema 15 before the read batch or report body is deserialized. They retain the ordinary 32 MiB envelope cap; the existing workspace-bearing 128 MiB whitelist is unchanged. Schemas 1–14 remain readable for the variants and semantics they introduced.

The request revision is `contig::batch_revision(&reads, &options)`: the exact ordered role IDs, complete molecular revisions, declared options and contig engine identity. Source evidence is bounded before hashing. The result's `batch_revision` must equal its response revision; report validators also bind the options and source/call evidence. Lightweight `validate_bindings_for` is appropriate after a worker reply; loading producer evidence never replaces a fresh assembly/reconstruction at a consensus-use boundary.

Recognized `read-assembly-v1` attachments inside older `SaveWorkspace` requests or `WorkspaceLoaded` replies require schema 15, even though those outer variants predate assembly. Native archive format 1 and plan schema 3 remain unchanged: assembly evidence is a separate retained producer record, not a new physical cloning operation. Candidate 5 and cloning verifier science-7 are unchanged.

Quality thresholds refer to their declared raw-score/explicit FASTQ model, not invented cross-format Phred probabilities. Branches, repeats, co-optimal placements, conflicts, unsupported gaps and work/output limits remain explicit. Model completion does not establish empirical uniqueness or independent physical read observations. [Scientific profile](CONTIGS.md).

## Methylation assessment (schema 16)

`AssessMethylation { molecule, options }` appends job 48; `Methylation(MethylationReport)` appends result 46. Both require schema 16 before their body is deserialized and retain the ordinary 32 MiB envelope cap. The existing workspace-bearing 128 MiB whitelist remains unchanged. All prior minimum-version guards, including the schema 15 contig attachment floor, remain in force.

The envelope revision is `methylation::batch_revision(molecule, options)`, binding the exact source, explicitly declared Dam/Dcm/CpG states, selected enzyme profiles, engine/profile identities and resource options. Source preflight precedes hashing. Response validation checks structural evidence and the batch binding; callers also check the expected immutable source/options using `validate_bindings_for`. Imported producer evidence requires fresh `replay` before use as reconstructed scientific evidence.

The assessment returns site/profile evidence only. It does not introduce a physical operation, digest product or verification pass. Candidate schema 5, workspace plan schema 3 and cloning science-7 remain unchanged. Unsupported profiles, ambiguous contexts, unknown/hemi/mixed source states and bounded incomplete searches retain explicit uncertainty. [Scientific profile](METHYLATION.md).

## Tracked generation and checked digest plans (schema 17)

`DeclareSourceModifications { template, assumptions }` and `CheckedDigest { template, enzymes, profile_revision, fragment_index }` append operation tags 11 and 12. Existing operation and report payload fields retain their layouts. New nested operations require worker schema 17 in verification/PDF jobs, retained-run requests, inline workspace graphs and workspace archive jobs/results. Science-8 verification/replay responses and new run records likewise require schema 17. Ordinary scientific messages retain 32 MiB bounds and the workspace whitelist retains 128 MiB bounds.

Candidate files write schema 6 and read schemas 1–5; tracked operations cannot be relabeled under an older candidate header. Native archives retain format 1 and write workspace plan schema 4, accepting historical plan schemas 1–3 subject to their operation floors. New retained runs use plan schema 4; old schema-3 run records remain historical evidence. New operations or science-8 assessments cannot hide in a schema-3 record. Existing report layouts remain readable with their original recorded engine identity; loading never changes science-7 or unrecorded legacy reports into current verification.

The generation declaration accepts only an explicit input source and declares Dam/Dcm/CpG generation assumptions in a closed three-mechanism model. A checked digest pins the separate positional sensitivity profile revision and consumes freshly reconstructed strand modification state. Serialized molecule summaries are producer evidence; they cannot restore tracked physical authority. All native/browser/MCP verification and protocol export rebuild the plan. Partial, unsupported, mixed or unknown effects cannot silently fall back to an unchecked digest.

## CRISPR computational designs (schema 18)

Jobs 49–52 append `DesignCrisprGuides { molecule, options }`, `SearchCrisprOffTargets { query, references, options }`, `DesignEditTemplate { source, target, options, guide }` and `CheckEditTemplate { source, target, options, guide, donor }`. Results 47–50 append `CrisprGuides(GuideReport)`, `CrisprOffTargets(OffTargetReport)`, `EditTemplate(EditTemplateReport)` and `EditTemplateCheck(EditTemplateCheck)`. Borrowed tag-prefix guards require schema 18 before any new variant body is decoded. All four retain the ordinary 32 MiB envelope cap; the existing workspace-bearing 128 MiB whitelist is unchanged.

Guide request revisions use `crispr::guide_revision(molecule, options)`. Off-target requests use `off_target_revision(query, references, options)`, binding the exact ordered reference roles/revisions/build labels/accessions, query, profile and all options. Edit-template requests use `edit_template::batch_revision(source, target, options, guide)`; independent checks use `check_batch_revision(..., donor)`, additionally binding the normalized concrete candidate sequence. Every response validates structural evidence and its envelope batch revision. Callers also use the matching lightweight current-engine/source/options binder; imported producer evidence never replaces fresh reconstruction/checking at use boundaries.

Source preflight precedes hashing. Native MCP reference batches use borrowed aggregate preflight before copying molecules. Off-target intended status is an exact role/revision/coordinate/orientation label, never a sequence-wide exclusion. Selected build labels cannot establish verified genome coverage. The editing target remains a noncausal requirement, and candidate checks reconstruct the expected edited source independently. These reports grant no physical editing operation or biological repair proof; candidate schema 6, workspace plan schema 4 and cloning science-8 remain unchanged. [CRISPR model](CRISPR.md); [editing-template model](EDIT_TEMPLATE.md).

## Streaming reference sessions (schema 19)

Jobs 53–66 append the stateful genome operations listed above. Results 51–59 append `GenomeImportStarted`, `GenomeImportProgress`, `GenomeFastaFinished`, `GenomeRegistered`, `GenomeIntervalGuides`, `GenomeSearchState`, `GenomeSearchPage`, `GenomeDisposed` and `GenomeIntervalRead`, using the DTOs in `genome_session`. Borrowed tag-prefix guards require schema 19 before decoding any new body. Every new job/result retains the ordinary 32 MiB envelope cap. Raw source chunks are at most 1 MiB; a native manifest is at most 16 MiB. The workspace-bearing 128 MiB whitelist and ordinary `Molecule` limits are unchanged.

`AnalysisSession` owns only parser/hash/current-block state, bounded manifest metadata, live search cursors and bounded hits. Hosts retain packed blocks in IndexedDB/native scratchfiles. FASTA push responses emit bounded `EncodedBlock` records for host staging. `FinishGenomeFasta` returns residual blocks and a `FinishedReference`, without semantic publication. The host reads every stored block back through `BeginGenomeNative`, ordered `PushGenomeNativeBlock` and `FinishGenomeNative`; only the private full validator can produce a fresh receipt and verified dataset registration. On reload, `RegisterGenome {manifest}` always creates unverified session state. Serialized receipts/progress cannot restore the private validator or search cursor.

Request revisions are:

| Jobs | Binding |
|---|---|
| 53 BeginGenomeFasta | `genome_session::fasta_revision(build, options, expected_source_bytes)` |
| 54/55/57/58/59 import push/finish/abort | Exact live `import_id` |
| 56 BeginGenomeNative | `native_revision(header, tail, file_length, manifest_bytes)` |
| 60 RegisterGenome | Exact `manifest.manifest_hash()` |
| 61/62/66 interval guides, search start, interval read | Exact `GenomeRef.manifest_hash` plus private handle/content lookup |
| 63/64 search advance/page | Exact live `job_id` |
| 65 CancelGenomeJob | Exact owned resource `id` |

The FASTA begin pins the declared raw source byte length. Chunks cannot exceed it, and finish rejects truncation before consuming the parser. Native framing pins its declared file length. Responses retain their matching transport binding and source/session hashes; `genome_session::validate_result` checks shape and binding without restoring authority. `ImportProgress.received` counts source bytes for FASTA, canonical bases for native readback.

A search advance must supply exactly the requested next owner plus halo blocks. Only the owner block enters ordered canonical-contig replay, avoiding double counting halos. Corrupt/skipped/replayed evidence cannot establish coverage. Whole-manifest absence evidence requires fresh private EOF validation and exhaustive model traversal; declared expected components and public build labels remain distinct. Interval guide/read reports attest local verified payloads only. Typed `genome-reference-link-v1` attachments require worker 19 inside older workspace variants but do not embed genome bytes or restore a session. Archive format 1, plan schema 4, candidate schema 6 and science-8 remain unchanged. [Large-reference model](GENOME_CRISPR.md).


## Desktop virtual digest (schema 20)

Jobs 67/68 append `InspectRestriction {source_role, molecule, enzymes, options}` and `MaterializeDigestFragment {source_role, molecule, enzymes, options, fragment_index}`. Result 60 is `RestrictionInspection`; result 61 is `DigestProduct`, using the DTOs in `restriction_inspection`. Both retain the ordinary 32 MiB envelope cap. Revision helpers preflight bounded rich source evidence before hashing. Inspection binding covers the exact stable role, source molecular revision, ordered enzyme names/recognition/cut geometry and all options; product binding additionally covers the zero-based source-cut-order fragment index. Response decoding validates shape and requires `batch_revision == response.revision`. These are cheap transport/current-input checks, not an independent scientific replay of imported reports.

Recognition scanning preserves exact recognition lengths and origin-split spans, possible matches and cuts outside the substrate. `scan_complete` and `products_determined` are separate. A completed recognition enumeration can still lack physical products; a limited scan cannot establish absence. Fragment summaries retain original top/bottom boundaries, lengths, ends and topology. Materialization freshly reconstructs the selected fragment, projects source annotations and removes obsolete original-record evidence. Producer recipes and duplex summaries remain descriptive; importing them cannot grant physical authority.

Job 69 appends `InstallDigestFragment {workspace, source, enzymes, options, fragment_index}`. `source` binds an exact document/event/revision through the existing history model. The request revision is `restriction_history::installation_revision`; it covers the workspace fingerprint, exact source event and selected product recipe. Result 62, `DigestInstalled`, returns a native archive and bound source/product identifiers. These two exact tags alone extend the workspace 128 MiB whitelist; they require schema 20 before a large body is decoded. Installation reruns the calculation against the current source and creates one derived document with the exact source event as parent. No serialized input product or new attachment kind is accepted. The receiver must verify current request/epoch bindings and validate the returned archive before installing it atomically. Candidate 6, plan 4, archive 1 and science-8 remain unchanged.


## Catalog-pinned digest (schema 21)

Jobs 70/71 append `InspectCatalogRestriction {selection, source_role, molecule, options}` and `MaterializeCatalogDigestFragment {selection, source_role, molecule, options, fragment_index}`. Job 72 appends `InstallCatalogDigestFragment {selection, workspace, source, options, fragment_index}`. Results 63–65 are `CatalogRestrictionInspection`, `CatalogDigestProduct` and `CatalogDigestInstalled`. Each is a wrapper from `catalog_digest`, retaining the original inspection/product/installation fields in its inner DTO. Jobs67–69/results60–62 and recipe1 remain readable without invented catalog identities.

A selection pins schema1, actual catalog edition/content hash and ordered exact enzyme IDs with full name/recognition/strand-cut geometry. Workers load the embedded sealed catalog themselves; clients cannot substitute arbitrary supplier records or promote a claimed hash. Revision helpers include this complete pin in a distinct domain around the existing source-role/geometry/options/index bindings. Before calculation, exact current catalog matching and bounded rich-source preflight are required. Response decoding checks wrapper shape plus its outer batch revision against the envelope. Imported historical reports remain producer evidence; lightweight validation does not replay biology or grant authority.

Only the exact installation request72/response65 tags extend the workspace128MiB whitelist, with a schema21 floor before large-body decode. Scientific requests70/71/results63/64 retain32MiB. Installation uses the existing fresh derivation/history pipeline to create one product with its exact source event as parent; it accepts no serialized product. `validate_archive_bindings_for` checks current request, decoded recipe and catalog bindings without running site matching on the main thread; it is used only after a fresh private worker result and current epoch guard. Native `validate_archive_for` additionally reconstructs the calculation. Source-only catalog comments cannot restore fresh chemistry. Candidate6, plan4, archive1, recipe1 and science-8 are unchanged.

## Reference recombination (schema 22)

Jobs 73–75 append `InspectRecombination {inputs, request}`, `CheckRecombinationTarget {inputs, request, target, policy}` and `SelectRecombinationComponent {inputs, request, selection, target, policy}`. Results 66–68 are `RecombinationInspection`, `RecombinationTargetCheck` and `RecombinationProduct`. All use the ordinary 32 MiB cap. Their request revisions are the adapter's `batch_revision`, `check_revision` and `checked_product_revision`, respectively. Source roles/order, full pinned profile snapshot, limits, annotation policy and independent target policy bind the relevant requests. Borrowed aggregate guards precede molecule copies, validation and hashing.

The inspection retains every reciprocal component and reviewed/unmapped source annotations. Target checking reconstructs sources afresh without selecting a component. Selected-product reconstruction requires the explicit candidate/component plus an independent target match. Lightweight DTO checks bind actual DNA, source identities, exact request/profile, annotation accounting and uncertainty; these checks do not confer positive scientific authority on imported reports. Fresh source-based reconstruction remains required.

`RecombineReference` is appended to the operation enum, requiring candidate 7, workspace plan 5, worker 22 and science-9. Historical candidate 1–6, workspace plan 1–4 and worker 1–21 layouts remain readable. New operations or engine claims cannot hide under older nested headers; the native archive container remains 1. Reference-only lineage stays chemically incomplete through descendants. Imported generated GenBank/EMBL markers can restrict interpretation but cannot restore a physical state. [Core model and bounds](RECOMBINATION_ADAPTER.md) · [Desktop workflow](RECOMBINATION_UI.md) · [MCP](RECOMBINATION_MCP.md).

## TOPO plan graphs (schema 25)

`Operation::Topo` appends operation tag 14, with explicit ordered source-role
bindings, a versioned generation/transfer request and a batch-bound component
selection. Existing Verify and ExportProtocolPdf jobs carry these plans; nested
TOPO operations require schema 25 even though their outer job tags are older.
Candidate schema 8 and workspace plan schema 6 prevent relabeling the operation
as an older plan. The physical graph is freshly reconstructed from source DNA,
not restored from a report. No unique paired contour means no invented reference
DNA. A contour preview retains restrictive chemistry boundaries; the initial plan
integration remains Incomplete without a physical acceptance/downstream model.

Verification reports append `step_reports` with complete source-backed TOPO
physical evidence. Worker schemas 3–24 decode their actual previous positional
report layout and gain an empty evidence list; schemas 1–2 additionally retain the
unrecorded engine marker. A historical report never gains a current graph or
current verification authority. `science-10` reports require schema 25.

SAM/BAM result tags 69–73 were appended in schema 24; these remain stable.
Read collections stay private to a worker session and decoded JSON or binary
receipts cannot reinstall parsed records or checked biological authority.

## Standalone TOPO inspection (schema 26)

Result tags 74–77 append fresh physical simulation, orientation inspection,
selected-component projection and independent target assessment. Request
revision helpers hash bounded original input identities, full recipe and limits,
ordered strand indexes or component/seed selection, and the actual target and
all acceptance-policy flags where applicable. Hashing a request performs no
scientific reconstruction. The worker compares that binding before running the
engine; decoding a report checks readonly shape and context, not biological
truth or experimental success. Reopening historical JSON cannot install graph
state. Browser source-event/epoch checks provide a separate identity boundary
for equal-content documents and changed histories.

## TOPO source annotations (schema 27)

Result tags78/79 append annotated projection and independent-target assessment.
The additive wrappers retain the unchanged physical Product plus separate
annotation evidence/reference. Hash domains include the annotation engine so an
old annotation-free projection cannot stand in for a newly requested view.
Reports append `annotation_reports`; explicit positional migration for schemas25
and26 retains their existing step graphs with an empty annotation-evidence list.
The current science11 verifier requires schema27. JSON defaults are for historical
representation only and cannot reconstruct proof or grant scientific authority.

## Final checked TOPO plans (schema 28)

Operation tag15 appends `CheckedTopo` with explicit source roles, preparation,
component/seed and all six acceptance flags. It is final-only and compares an
independent bound target via fresh physical reconstruction. Older job tags
Verify/PDF cannot hide the operation beneath a pre28 envelope. Candidate files
require9, workspace plans require7, and science12 reports require28. No new
report field or standalone worker job is introduced, preserving the schema27
report layout. Explicit reverse/rotation policies are also used by native run
capture; acceptance never promotes a reference view into an intact duplex.

## Variant sessions (schema29)

Results80–85 append completed-source summaries, bounded normalized pages, fresh
actual-reference checks, fresh selected-allele previews, exports and disposal.
Source hashes and complete request options bind every result. A private parsed
Collection is required; JSON summaries and decoded reports cannot reconstitute
one. Stateless routing rejects session-required jobs. Preview DNA retains its
chemistry-unmapped design guard and is not a confirmed genotype or prepared
substrate. Explicit reference-event identity adds a separate browser lifecycle
boundary, including equal-content documents and workspace replacement.

## GTF annotation selector (schema 30)

The append-only `AnnotationFormat::Gtf` selector is 2; GFF3 remains 0 and BED remains 1. Existing annotation jobs 24–26 and import report tag retain their layouts. GTF import/export, applying a GTF report and receiving a GTF import report require schema 30; rewriting their envelope as schema 29 or earlier is rejected. GFF3/BED retain their schema-7 contract. Retained source text remains unverified annotation evidence; applying it freshly reconstructs the report against the chosen reference.

## Prepared overlap assembly (schema 31)

Job 94, `PreviewOverlapRecipe`, accepts immutable original `PlanInput` roles, an
independently supplied target and an explicit ordered recipe. Result 86,
`OverlapRecipePreview`, retains the compiled ordinary plan, freshly reconstructed
preparations, verifier output and exact/advisory overlap evidence. PCR, digestion
and optional end repair run in the worker; a document excerpt is not treated as
physical preparation. Whole-source use requires a linear source.

The context binds source IDs/revisions, target revision and every recipe option.
The new outer selectors require schema 31 before payload decoding. Existing
candidate (9), workspace-plan (7) and scientific-engine (science-12) versions stay
unchanged because the compiler emits existing operations. Imported preview data
cannot recreate private chemistry or bypass a fresh plan check. Ordinary 32 MiB
message limits apply; no oversized workspace-message exemption is added.

## Prepared ordered ligation (schema 32)

Job 95, `PreviewLigationRecipe`, carries actual original `PlanInput` roles, an
independent target and ordered preparation/treatment/orientation/nick policy.
Result 87, `LigationRecipePreview`, retains the compiled existing-operation plan,
fresh verifier results, oriented preparation-end summaries, junctions and actual
product nick state. The physical summaries come from private material state in
the same verifier execution, not from deserialized product metadata. Unknown
terminal chemistry remains unknown; nick acceptance does not invent a bond.

Context binds all source/target revisions and recipe options. Selectors require
schema 32 before payload decoding. Existing overlap schema-31 DTOs are unchanged;
candidate 9, workspace-plan 7 and science-12 remain unchanged because both recipe
compilers emit existing scientific operations. Ordinary 32 MiB limits apply.

## Assembly junction inspection (schema 33)

Job 96, `InspectAssemblyJunctions`, accepts the actual typed plan, independently
bound target and explicit flank/genetic-code/frame/inspection limits. Result 88,
`AssemblyJunctionInspection`, carries fresh verification, reconstructed product
identity, bounded actual DNA windows, shared-homology intervals and requested
window-relative peptide/CDS evidence. Final reconstruction must have succeeded
before a product can be inspected; an earlier preparation cannot substitute.
A completed product that disagrees with its independent target remains visible
alongside the invalid verification outcome.

The context binds the whole original plan, actual target and options. Both new
selectors require schema 33 before decoding their scientific payloads. Existing
plan/candidate/report layouts remain unchanged; native candidate 9, workspace-plan
7 and science-12 still apply. The ordinary 32 MiB message budget is retained.
The report is readonly evidence and cannot restore chemistry or verification
from imported metadata. Frames/translation do not establish expression, reading
frame intent or biological function.

## Explicit original-source termini (schema 34)

`Operation::DeclareSourceTermini` is appended at operation tag 16. It declares
left/right `Terminus` hypotheses on an eligible original intact linear DNA input,
without altering source DNA or inferring chemistry from a file. `Preparation::
DeclaredDuplex` is appended at preparation tag 3 and compiles to that operation.
Existing enum tags and serialized struct layouts are unchanged.

New semantics require worker 34, candidate 10, workspace-plan 8 and science-13.
Checks cover direct Verify/PDF requests, the existing preparation preview jobs,
nested candidate graphs, inline/native workspace jobs and retained-run records.
A newer declaration or engine report cannot be relabeled under an older header.
Historical candidates 1–9, workspace plans 1–7 and worker messages 1–33 remain
readable within their original semantics; loaded reports never restore chemistry.
The outer native archive/document formats remain version 1.

Original physical/reference-only markers and ssDNA/circular sources remain
restrictive. Declaring ends on a step is unsupported, unknown groups stay unknown,
and a declaration neither repairs material nor creates a methylation-generation
claim. HTML/PDF retain the explicit operator hypothesis. [Scientific contract](DECLARED_SOURCE_TERMINI.md).

## Source-derived overlap primer suggestions (schema 35)

Append-only job 97 (`DesignOverlapPrimers`) and result 89
(`OverlapPrimersDesigned`) bind original source IDs/revisions, explicit source
intervals/orientations, primer/profile settings and the separate target in a
context revision. Both tags reject pre-35 headers before deserializing the new
payload. The ordinary 32 MiB cap applies. This is an advisory designer, not a
verification report; applying a suggestion produces ordinary exact-PCR recipe
rows for the existing fresh compiler/verifier and protocol exports. No new
candidate, workspace-plan or scientific-operation layout is introduced.

## Frozen assembly design intent (schema 36)

Append-only job 98 (`BuildAssemblyIntent`) and result 90 (`AssemblyIntent`) build
an expected reference from explicit source regions, order and orientation. The
request contains actual original inputs and a region request, with no complete
target or reaction operation. Context binds source IDs/revisions and the request;
changed regions/topology/source revisions fail stale requests. New tags reject
pre-36 headers before their payload is decoded, under the ordinary 32 MiB cap.
The response describes declared intent and contribution coordinates, not a
scientific verification result. Its reference-only target is separately passed
to existing preparation/planning and fresh verification. Candidate 10, plan 8
and science-13 layouts are unchanged.

## Catalog-pinned digest preparation (schema 37)

The append-only operation `CatalogDigest` is tag 17 and shared preparation
`CatalogDigest` is tag 4. They carry an exact catalog selection and fragment index;
preparation can append explicit end repair. Fresh reconstruction validates the
selection against the embedded catalog and uses the same tracked digest engine
as existing restriction operations, without converting catalog IDs to starter-only
names. Complete cleavage remains an explicit assumption.

Existing overlap/ligation recipe jobs require schema 37 when their rows contain
this preparation. Plans with the new operation require worker 37, candidate 11
and workspace-plan 9. Science-14 reports also require worker 37; nested candidate,
run and archive semantics are checked even when wrapped in an older message.
Historical operations and reports retain their previous minimum versions.

## Direct cut selection and reviewed outputs (schema 38)

`SelectDigestFragment` (job99/result91) binds the original source role/revision,
full catalog selection, inspection options, expected fresh catalog batch and two
explicit boundaries. It reruns inspection and materialization, accepts one actual
linear fragment with matching top/bottom intervals, and refuses paths split by
additional cuts. Circular complement exchanges endpoints; identical circular cut
instances select one full turn. Linear source ends are explicit boundaries.

`ExportPlanOligos` (job100/result92) binds the exact plan, independent target,
deduplication option and CSV/TSV delimiter. It extracts declared order data without
rerunning PCR or assigning a positive verification status. The typed result
includes the table, every use's domains/tails/policies, unavailable diagnostics and
spreadsheet-safe text. Complete=false is not an all-required order sheet.

`ExportProtocolHtml` (job101/result93) uses the same bound plan/target as the
existing PDF exporter, with fresh scientific reconstruction in the worker. Its
response is escaped HTML text under the ordinary message cap. New request/result
tags refuse pre-38 headers before payload decoding. These additions do not change
candidate11, plan9 or science-14 operation layouts.

## Atomic cloning outputs (schema 39)

`FinalizeCloningOutputs` (job 102/result 94) carries a complete workspace, exact
source/target event bindings, a plan and explicit output options. It reconstructs
the run independently, preserves every preparation event, optionally saves a
complete primer collection, and records presentation-only hidden document IDs.
Its context revision includes the workspace, bindings and options. Receipts validate
that revision and the returned workspace relationships before UI installation.
These exact workspace-bearing tags have the bounded 128 MiB limit and reject
pre-39 headers before payload decoding; ordinary jobs retain 32 MiB.
Candidate 11, plan 9 and science-14 are unchanged. See [outputs](CLONING_OUTPUTS.md).

## PCR-added Type IIS sites (schema 40)

`DesignGoldenGatePrimers` (job103/result95) binds original immutable inputs, an
independent target and explicit source intervals/orientations/domain lengths/right
fusion DNA, enzyme, clamp and spacer. The report contains source-bound segmented
oligos and `PcrCatalogDigest` preparation requests. It does not install a plan or
verify the intended target. New tags reject pre-40 headers before payload decoding;
request and response use the ordinary 32 MiB cap.

Shared preparation `PcrCatalogDigest` is append-only tag5. Existing overlap/ligation
recipe request tags require schema40 when their nested rows use it; old headers
cannot hide the new preparation. Compilation emits the existing PCR and pinned
catalog digest operations, plus optional explicit repair, so candidate11, plan9 and
science-14 do not change. `PreparedRow.preparation_step` is the first PCR step,
while `prepared_source` points to the final digest/repair result. No existing row
fields or historical layouts change. See [planner](GOLDEN_GATE_PRIMER_DESIGN.md).

## Explicit mutation intent (schema 41)

`BuildMutationIntent` (job104/result96) binds an original source molecule and an
explicit named replacement interval/inserted DNA request. The report contains the
unverified intended reference, deleted source bases, actual source revision and
explicit circular-origin shift. It is independent of all preparation simulators.
New tags reject pre-41 headers before payload decoding and use ordinary32MiB limits.
The target carries restrictive reference-only metadata, including source/request/
sequence receipts. Readonly report validation does not establish prepared material.

Candidate11, plan9, science-14 and native archive1 remain unchanged. The UI retains
exact source events outside the scientific DTO, checks them before opening, and
records computational parentage and `mutation-intent-v1` evidence. Inverse-PCR
planning separately reconstructs against this intended target, including the valid
circular terminal insertion caret. See [intent scope](MUTATION_INTENT.md).

## Construct-owned read evidence (schema 42)

`CaptureReadCheck` (job105/result97) carries the full workspace, exact construct
reference/read event bindings and native validation options. It resolves retained
events, runs fresh validation and returns a private workspace plus immutable saved
check. The installation context binds the full incoming workspace; the stored
check context binds events/options/engine, permitting content-addressed retry dedupe.
`ReplayReadCheck` (job106/result98) resolves the saved original events and returns
both historical record and newly reconstructed current-engine evidence.

Both new jobs and the capture workspace result use the explicit128MiB workspace
whitelist. The ordinary replay result retains32MiB. New tags refuse pre-42 headers
before decoding; nested workspace semantics remain checked. Loading a record never
restores scientific authority, and partial coverage cannot become whole-construct
confirmation. Candidate11/plan9/science-14/archive1 are unchanged.
[Read-check contract](READ_CHECKS.md).

## Explicit review decisions (schema 43)

`PreviewReadReview` (job107/result99) and `CommitReadReview` (job108/result100)
carry the native workspace and an explicit source-check-bound decision request.
Each operation freshly replays original reference/read events. Nominated bases
include an expected concrete allele plus read/call identity; insertions include
expected read, boundary and full sequence. Changed evidence meaning cannot silently
retarget an index-only choice. Retain/unresolved choices do not change DNA.

The two workspace-bearing requests and committed workspace result use the explicit
128MiB whitelist. The ordinary preview result retains32MiB. New tags reject pre-43
headers before decoding, and nested workspace semantics remain checked. Result
contexts bind the complete incoming workspace/request. Saved records remain
historical; the optional target is a reference proposal with exact source parents,
not a restored current read verdict. Candidate11/plan9/science-14/archive1 are
unchanged. [Review behavior](READ_REVIEW.md).

## Junction-edited intent and oligos (schema 44)

`BuildEditedAssemblyIntent` (job109/result101) binds immutable original inputs,
the original region/order/orientation request and explicit junction trims/inserts.
It returns the independent edited reference and adjusted original-source coordinates;
no PCR/overlap product supplies expected DNA. Output contribution spans include each
core's declared outgoing synthetic suffix. `DesignEditedAssemblyPrimers`
(job110/result102) additionally binds overlap/profile/primer controls and proposes
explicit original-source PCR rows with the full insert-plus-homology tails analyzed
before ranking. It remains advisory until fresh recipe reconstruction.

Both request/result pairs use ordinary32MiB limits and reject pre-44 headers before
decoding. The old pure region-intent DTO/report layout is unchanged; an edited target
is a distinct typed report, not a rewritten old report. Candidate11, plan9,
science-14 and archive1 remain unchanged. [Edit geometry](ASSEMBLY_JUNCTION_EDIT.md).


## Selective junction proposals (schema 45)

Append-only job111 `DesignSelectiveAssemblyPrimers` returns result103
`SelectiveAssemblyPrimers`. The request binds original sources, base intent,
junction edits, ordered Retain/DesignPcr rows, explicit per-junction overlaps and
primer/profile options. Retain preserves the full preparation DTO; only selected
rows receive new original-source oligos. The independent edited target is built
before fresh complete-recipe verification, which remains a separate operation.
Both directions require schema45 before payload decoding and validate the request
context. Ordinary32MiB limits apply; no workspace-size exception is introduced.
Candidate11, plan9, science14 and archive1 remain unchanged.


## Golden Gate junction proposals (schema 46)

Append-only job112 `DesignGoldenGateJunctionPrimers` returns result104
`GoldenGateJunctionPrimers`. The request binds original source events via exact
input revisions, base declarations, junction trims/inserts, retained or explicitly
selected PCR rows, enzyme/clamp/spacer and terminal/nick policy. The target is the
independent edited source-region intent; it is never used as an amplification
template. New own-5′ primer segments distinguish fusion and inserted bases.
Selected originals undergo fresh PCR, actual pinned Type IIS cut and oriented
end checks before returning an advisory complete ordered-ligation recipe.
Retained rows stay exact and can still fail a subsequent complete-recipe check.
Both directions require schema46 before payload decoding; ordinary32MiB caps apply.
Candidate11, plan9, science14 and archive1 are unchanged. Old unedited Golden Gate
DTOs, tags and geometry remain unchanged.


## Prepared Golden Gate pool inspection and selection (schema 47)

Append-only job113 `InspectPreparedGoldenGatePool` returns result105
`PreparedGoldenGatePool`. It binds original inputs, an independent target and
explicit preparation recipe/pool options. Preparation and raw physical substrates
come from one fresh native pass; missing required material withholds the complete
pool rather than silently analyzing an available subset. The existing finite
Golden Gate engine redigests those selected prepared roles, reports backgrounds,
recutting candidates, caps and actual target correspondence. It does not cover
unrepresented discarded fragments or all material in an original one-pot reaction.

Job114 `SelectPreparedGoldenGateCandidate` carries a boxed bound report and exact
parent context/index/directed+physical keys, returning result106
`PreparedGoldenGateCandidate`. Selection preserves original inputs and preparation
steps and replaces only final ordered ligation with the existing GoldenGate
operation. It returns an unverified plan, never trusted physical state or product
DNA. The standard verifier must freshly reconstruct that plan. Request/result
bindings include the derived plan, selector and parent context; both directions
require schema47 before payload decoding. Ordinary32MiB limits apply throughout.
Candidate11, plan9, science14 and archive1 remain unchanged.


## Selected Golden Gate junction windows (schema 48)

Existing `InspectAssemblyJunctions` job96/result88 now supports final GoldenGate
operations using a private enumeration receipt captured in the same fresh verifier
pass. The appended nested method selector3 is `GoldenGate`, with actual path/source
mapping and finite-search metadata. Earlier method selectors0–2 and their payloads
are unchanged. Response decoding rejects this new method under pre48 headers.
The existing request payload remains valid at its original version: it already
carried GoldenGate plans, previously reported as unsupported junction geometry.
No new top-level job/result, operation or large-message exception is introduced.
Ordinary32MiB caps and candidate11/plan9/science14/archive1 remain unchanged.


## Explicit reaction substrates (schema 49)

Append-only inspect115/result107 binds originals, independent target and explicit
WholeSource/Prepared reaction roles. It captures the actual final enumeration in
one verifier pass, including unavailable or recuttable probe index0, without
publishing a candidate. Select116/result108 binds an exact context/index/key pair
into an unverified native plan. It uses the existing generic key selection/report
shape while keeping a separate whole-reaction context and result selector.
Pre49 envelopes reject both new directions before payload decode. Ordinary32MiB
limits apply; old prepared-pool DTOs, candidate11, plan9, science14 and archive1
remain unchanged.


## Selected PCR reference outputs (schema 50)

Job117 `PreparePcrOutputs { request }` binds the full original source, named
oligos, selected geometry index/preview revision and explicit output choices.
Result109 `PcrOutputs` retains the freshly reconstructed geometry and optional
named reference/oligo collection. Both new selectors require schema50 before
payload decoding and use the ordinary 32 MiB envelope. The host separately binds
source events and the live draft before installing outputs. This adds no physical
operation or independent-target verification. [Output contract](PCR_OUTPUTS.md).

## Selected conditional mismatch-PCR outputs (schema 51)

Job118 `PrepareMismatchPcrOutputs { request }` binds the original source, named
pair, actual mismatch/budget options, selected original preview revision and
output choices. Result110 `MismatchPcrOutputs` retains the freshly selected native
geometry, optional reference/collection and full search summary. Both selectors
require schema51 before payload decoding; ordinary 32 MiB limits apply. No new
scientific operation or candidate format is introduced. Hosts retain exact source
history and live draft checks. [Output contract](PCR_OUTPUTS.md#conditional-mismatch-pcr-outputs).

## Event-bound pairwise comparison (schema 52)

Job119 `CompareSequences { reference, query, options }` binds both full source
revisions and either stored-orientation affine options or explicit whole-molecule
origin/strand budgets. Result111 `ReviewedComparison` retains the native result,
source lengths, context and interpretation version. Both selectors require
schema52 before new payload decoding and use the ordinary 32 MiB envelope.
Existing Align job/result layouts remain supported. Exact history events and
current/historical state are host/archive concerns; neither checksum nor saved
report proves native optimality. [Comparison contract](PAIRWISE_REVIEW.md).

## Annotated reference transfer (schema 53)

Jobs120 `CopyReferenceFragment { request }` and121
`PasteReferenceFragment { request }` return112 `ReferenceFragment` and113
`ReferencePaste`. Copy context binds the actual source and ordered selection;
paste binds actual original source/destination, copied context and destination
interval. The native paste rebuilds the fragment instead of trusting serialized
clipboard letters or annotations. Both new selectors require schema53 before
payload decoding and retain ordinary 32 MiB envelopes. No new physical operation
or workspace-size permission is introduced. [Transfer contract](REFERENCE_TRANSFER.md).


## Source-bound translated references (schema 54)

Append-only jobs122 `TranslateReference { request }` and123
`PrepareTranslatedReference { request, name }` return114 `TranslatedReference`
and115 `TranslatedOutputs`. Requests bind the complete original source plus an
explicit compound region/orientation/frame/code policy or an exact CDS feature
index. Inspection returns the native translation report, protein letters,
effective location and any CDS/qualifier diagnostics. Output preparation freshly
reconstructs the same translation and returns a named linear protein reference
and, when supported, an explicit source annotation proposal. Output names are
part of their separate context revision; inspection context remains separate.

Both directions reject pre54 envelopes before decoding their payloads; schema53
reference-transfer selector gates remain53. These are ordinary scientific jobs:
bounded serialization and response/request decoding retain the 32 MiB envelope,
without the workspace exception. Native source preflight and context checks run
before translation; response DTO validation checks its expected context and
bounded native content. The host still binds exact source events and its current
draft, and must explicitly apply outputs/annotation proposals. Imported protein
letters or a serialized report do not authorize source changes. No physical
operation, candidate version, genetic-code expansion or alignment executable is
introduced; candidate11, plan9, science14 and native archive1 remain unchanged.

## Catalog-pinned predicted gels (schema 55)

Job124 `PredictPinnedGel { molecule, lanes, calibration }` returns116 `PinnedGel`.
Pinned lanes carry an exact source event and full catalog selection/inspection
options, then freshly inspect actual source calls. Existing `PredictGel`/`Gel`
layouts remain readable; appending `Preparation::CatalogDigest` does not alter
fragment or report fields. An old `PredictGel` operation cannot transport a pinned
recipe, and pre55 prefixes cannot relabel the new job/result. Ordinary 32 MiB
message caps apply; no workspace-size allowance is added.

The host resolves historical source ownership from its retained workspace without
moving document heads. Native prediction validates actual molecular revision and
catalog geometry, requires complete/determined cleavage, and retains unplaced
circular conformations and explicit lane errors. Saved reports are documentary
evidence, never substitutes for fresh source reconstruction.

## Computed comparison presentation (schema56)

Job125 `ExportPairwiseComparison` returns117 `ExportedPairwiseComparison` for
aligned FASTA/Clustal or a complete HTML/PDF review. Job126
`SearchPairwiseComparison` returns118 `SearchedPairwiseComparison`; literal and
compatible-IUPAC searches retain gap-column windows and original source mappings,
explicit work/hit limits and incomplete status. Job127 `PreparePairwiseContext`
returns119 `PreparedPairwiseContext`, caching actual-source annotation and full-CDS
mapping for paged projection without restarting frames or doing translation on
the GUI thread. [Contracts and limits](PAIRWISE_PRESENTATION.md).

Each request binds all relevant inputs to its exact revision and stays within the
ordinary32 MiB cap. Old55 envelopes reject the new selectors before payload decode,
including short-prefix relabeling. Result guards bound byte/coordinate/cache data;
the host separately checks exact events, request/epoch/generation and current form.
These are views of an existing comparison, not new alignment or verification
authority. Existing job/result layouts and native/candidate/plan/science versions
are unchanged.


## Atomic combinatorial release and archives (schema58)

Begin131 owns a private workspace and selected historical declarations. Advance132
reconstructs/finalizes one row; Inspect133 exposes actual staged product and target
molecules; Finish134 returns all outputs only after every row succeeds. Any failed
row poisons the session. Responses123–126 carry started/advanced/inspected/finished
DTOs; serialized progress or manifests cannot restore a private session.

Export135/Archive127 uses exact released native events, an explicit format/subset
request and a bounded ZIP. HTML/PDF requests perform fresh current-engine checks;
sequence/primer-only exchange remains labelled historical. Begin/Finish/Export
workspace requests and Finish/Archive responses have explicit128 MiB envelope
allowances. Advance/inspection/progress retain32 MiB. Pre58 labels reject new
variable-length selectors before payload decoding. Existing candidate11, plan9,
science14 and native archive1 remain unchanged.
[Release and archive contracts](COMBINATORIAL_RELEASE.md).


## Exact Golden Gate candidate identity (schema59)

No outer job/result tags are added. The appended operation18,
`SelectedGoldenGate`, requires worker59/candidate12/plan10 and carries explicit
candidate index plus directed and physical sequence keys. One fresh native pool
enumeration checks those keys before accepting the product; the independent target
is checked separately. Older index-only Golden Gate operations retain their
original selector and readable historical layouts.

Nested `GoldenGatePool` method declarations in batch requests, including per-row
overrides, require59. Older ordinary batch payloads retain their57/58 contract.
New scientific reports identify science15 and cannot be relabelled58; workspace
candidates, recorded runs and batch/release attachments also preserve their nested
minimum versions. Saved candidate intent never restores a private session or
fresh authority. Envelope limits are unchanged. [Reaction matrix contract](COMBINATORIAL_POOLS.md).


## Published ligation profiles and physical competitors (schema62)

Jobs143 `ListLigationProfiles`,144 `ScorePublishedLigation`,145
`RankGoldenGateFusions` and146 `ScoreLigationCompetition` append results135–138.
New tags reject pre62 envelopes before decoding. All use the ordinary32 MiB cap;
no large-workspace allowance is added. Profile scoring binds exact profile/model,
intended junctions and separately declared literal physical competitors. Ranking
binds all original sources and the explicit conditional formulation hypothesis.
Catalog/revision/result guards are cheap and do not decode XLSX on the main thread.

Free linear outer ends enter denominators only, in their actual5′→3′ orientation;
absent complements and intended joins are not invented. The original job138 /
result130 supplied-count layout and model remain unchanged. Candidate12/plan10 /
science15 are unchanged. [Scientific and UI contract](PUBLISHED_LIGATION_PROFILES.md).


## Downstream TOPO nick sealing (schema63)

Job147 `SealTopoNicks` returns139 `TopoNicksSealed`. It binds original source
molecules, preparation/transfers and explicit downstream ligase/nick choices.
Pre63 envelopes reject new selectors before payload decoding. Ordinary32 MiB
message limits apply. Decoded evidence is structurally bound; scientific authority
requires fresh reconstruction, never a caller-supplied graph or positive receipt.
Existing operation/plan/candidate/science layouts remain unchanged.
[Scientific contract](TOPO_NICK_SEALING.md).


## Fresh processed TOPO contours (schema64)

Job148 `ProjectProcessedTopo` returns140 `ProcessedTopoProduct`, binding original
sources/preparation, explicit nick treatment, component/source seed and proof
options. Pre64 envelopes reject these selectors before payload decode. Ordinary
message caps apply; no workspace-size allowance is added. Full original/treated
graphs and conservative annotation/state evidence remain retained; imported result
validation cannot restore a fresh physical material. Existing old DTO/operation /
plan/candidate/science layouts remain unchanged.
[Contract](TOPO_PROCESSED_PRODUCT.md).


## Processed cloning plans and historical graph download (schema65)

Operation19 ProcessedTopo requires worker65/candidate13/plan11. Science16 reports
append processed_step_reports without changing old StepReport. Explicit legacy
<=64 decoders preserve all nine direct/nested VerificationReport paths; existing
1–26 bare-report converters run first. Missing historical evidence is never
invented. Nested plans and current science/evidence cannot be relabelled old.

Job149 LoadProcessedRunEvidence returns141 HistoricalProcessedRunEvidence,
containing exact historical JSON bytes and bounded run/engine/step metadata. Only
its workspace-bearing request joins the128 MiB whitelist; the result retains the
ordinary32 MiB cap with explicit export-workspace fallback at the boundary.
Historical graphs cannot restore fresh material or verification authority.
[Retained-plan contract](TOPO_PROCESSED_PLAN.md).


## Ordered physical processing (schema66)

Operation20 `TopoActions` requires worker66/candidate14/plan12/science17.
Jobs150 `ProcessTopoActions` and151 `InspectTopoActions` return142
`TopoActionsProduct` and143 `TopoActionsInspected`. Full typed evidence uses the
appended generic `physical_processing_reports` vector. Explicit nine-path
positional migrations through65 preserve old nonempty processed reports. Job149
also downloads a bounded exact-evidence bundle for generic/mixed retained runs.
[Chemistry, UI, MCP, historical binding and limits](TOPO_PROCESSING_ACTIONS.md).


## Upstream PCR template composition (schema67)

`TopoActions` with a `Source::Step` requires worker67/candidate15/plan13 and a
fresh native PCR-only template check. Science18 reports require67. No DTO layout,
operation tag, report field or worker selector changes. Original-input operation20
retains its66/14/12 floor and old reports remain documentary data.
[Native authority, synthesis and exact history](TOPO_UPSTREAM_PCR.md).
