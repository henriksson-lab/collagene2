# Binary worker protocol

The web worker and native adapters share Rust DTOs in `core::protocol` and the
same deterministic scientific engines. Current schemas are worker **123**,
candidate **31** and native workspace plan **29**; the verifier is **science-34**.
Only current internal encodings are accepted. Unsupported development versions
have no decoder, migration or capability ladder.

## Envelope and validation

Each frame has a 42-byte header: four-byte magic (`CGQ\0` request, `CGR\0`
response, `CGP\0` candidate), little-endian u16 schema, little-endian u32 payload
length, and 32-byte BLAKE3 payload checksum. The remaining payload is Postcard.
A worker request carries an ID, exact caller revision/context and a typed job.
Its response retains the ID and revision with a typed result. Consumers still
validate their captured source/event/workspace state before accepting a response.

The decoder checks namespace, current version, exact length, checksum and complete
payload consumption. Ordinary messages have a **32 MiB** cap; only explicitly
listed compound jobs/replies may use **128 MiB**. A borrowed prefix probe
checks the tag and bounded revision before admitting a large body. The current request selector list is shared by stateless and session admission;
reply admission has its own current selector list. These bounds are authoritative;
adding a job does not grant it the larger budget automatically.

Typed validators preserve scientific model bounds, exact source revisions, graph
dependencies and receipt correspondence. Nested material origins share ancestry,
byte, base and step budgets and are charged before decoding. A checksum establishes
byte integrity; it grants neither permission nor physical-material authority.
Current external biological import formats are separate from internal schemas.

## Scientific boundary

Plans declare ordered operations against an independently supplied intended target.
The native verifier reconstructs those operations using algorithms, never an AI
assertion. One `Ligate` operation retains its intended path and explicitly bounded
alternative products. A successful transport, saved assessment or documentary
material summary does not substitute for fresh native reconstruction.

Saved-run review resolves original source events and the retained independent
target without changing document heads. New protocol export uses fresh verification;
recorded findings remain documentary. SVG map export preserves the exact source
revision and supplied coordinate evidence, with explicit presentation limits.

## Current job tags

These are indexes in the current Postcard enum. The current schema must change
when an internal enum layout changes.

| Tag | Job |
| ---: | --- |
| 0 | `Search` |
| 1 | `Restriction` |
| 2 | `Pcr` |
| 3 | `Verify` |
| 4 | `DuplexDigest` |
| 5 | `Translate` |
| 6 | `TranslateCds` |
| 7 | `FindOrfs` |
| 8 | `Align` |
| 9 | `PcrDocuments` |
| 10 | `AnalyzePrimer` |
| 11 | `DesignPrimers` |
| 12 | `ExportProtocolPdf` |
| 13 | `GoldenGate` |
| 14 | `MismatchPcr` |
| 15 | `InversePcr` |
| 16 | `PredictGel` |
| 17 | `DiagnosticDigests` |
| 18 | `MismatchPcrPreview` |
| 19 | `SaveWorkspace` |
| 20 | `LoadWorkspace` |
| 21 | `ValidateTraces` |
| 22 | `ScanFeatureLibrary` |
| 23 | `ApplyFeatureLibrary` |
| 24 | `ImportAnnotations` |
| 25 | `ApplyAnnotations` |
| 26 | `ExportAnnotations` |
| 27 | `RepairPreview` |
| 28 | `AnnealPreview` |
| 29 | `AnnealedDocument` |
| 30 | `ViewTracks` |
| 31 | `SearchDocument` |
| 32 | `PrepareTracks` |
| 33 | `CachedViewTracks` |
| 34 | `ImportAlignment` |
| 35 | `ExportAlignment` |
| 36 | `ExtractAlignmentRows` |
| 37 | `SaveAlignment` |
| 38 | `LoadAlignment` |
| 39 | `ProteinProperties` |
| 40 | `TailPreview` |
| 41 | `CheckOverlapProfile` |
| 42 | `RecordCloningRun` |
| 43 | `ReplayCloningRun` |
| 44 | `DesignCodons` |
| 45 | `CheckCodonDesign` |
| 46 | `DotPlot` |
| 47 | `AssembleContigs` |
| 48 | `AssessMethylation` |
| 49 | `DesignCrisprGuides` |
| 50 | `SearchCrisprOffTargets` |
| 51 | `DesignEditTemplate` |
| 52 | `CheckEditTemplate` |
| 53 | `BeginGenomeFasta` |
| 54 | `PushGenomeFasta` |
| 55 | `FinishGenomeFasta` |
| 56 | `BeginGenomeNative` |
| 57 | `PushGenomeNativeBlock` |
| 58 | `FinishGenomeNative` |
| 59 | `AbortGenomeImport` |
| 60 | `RegisterGenome` |
| 61 | `GenomeIntervalGuides` |
| 62 | `StartGenomeSearch` |
| 63 | `AdvanceGenomeSearch` |
| 64 | `PageGenomeSearch` |
| 65 | `CancelGenomeJob` |
| 66 | `ReadGenomeInterval` |
| 67 | `InspectRestriction` |
| 68 | `MaterializeDigestFragment` |
| 69 | `InstallDigestFragment` |
| 70 | `InspectCatalogRestriction` |
| 71 | `MaterializeCatalogDigestFragment` |
| 72 | `InstallCatalogDigestFragment` |
| 73 | `InspectRecombination` |
| 74 | `CheckRecombinationTarget` |
| 75 | `SelectRecombinationComponent` |
| 76 | `ScanRestriction` |
| 77 | `PrepareReadAlignments` |
| 78 | `ReadAlignmentPage` |
| 79 | `CheckReadAlignment` |
| 80 | `ExportReadAlignments` |
| 81 | `DisposeReadAlignments` |
| 82 | `SimulateTopo` |
| 83 | `InspectTopoOrientations` |
| 84 | `ProjectTopoProduct` |
| 85 | `CheckTopoTarget` |
| 86 | `PrepareVariants` |
| 87 | `VariantPage` |
| 88 | `CheckVariantReference` |
| 89 | `PreviewVariantAllele` |
| 90 | `ExportVariants` |
| 91 | `DisposeVariants` |
| 92 | `PreviewOverlapRecipe` |
| 93 | `PreviewLigationRecipe` |
| 94 | `InspectAssemblyJunctions` |
| 95 | `DesignOverlapPrimers` |
| 96 | `BuildAssemblyIntent` |
| 97 | `SelectDigestFragment` |
| 98 | `ExportPlanOligos` |
| 99 | `ExportProtocolHtml` |
| 100 | `FinalizeCloningOutputs` |
| 101 | `DesignGoldenGatePrimers` |
| 102 | `BuildMutationIntent` |
| 103 | `CaptureReadCheck` |
| 104 | `ReplayReadCheck` |
| 105 | `PreviewReadReview` |
| 106 | `CommitReadReview` |
| 107 | `BuildEditedAssemblyIntent` |
| 108 | `DesignEditedAssemblyPrimers` |
| 109 | `DesignSelectiveAssemblyPrimers` |
| 110 | `DesignGoldenGateJunctionPrimers` |
| 111 | `InspectPreparedGoldenGatePool` |
| 112 | `SelectPreparedGoldenGateCandidate` |
| 113 | `InspectGoldenGateReaction` |
| 114 | `SelectGoldenGateReaction` |
| 115 | `PreparePcrOutputs` |
| 116 | `PrepareMismatchPcrOutputs` |
| 117 | `CompareSequences` |
| 118 | `CopyReferenceFragment` |
| 119 | `PasteReferenceFragment` |
| 120 | `TranslateReference` |
| 121 | `PrepareTranslatedReference` |
| 122 | `PredictPinnedGel` |
| 123 | `ExportPairwiseComparison` |
| 124 | `SearchPairwiseComparison` |
| 125 | `PreparePairwiseContext` |
| 126 | `BeginCombinatorialBatch` |
| 127 | `AdvanceCombinatorialBatch` |
| 128 | `CaptureCombinatorialBatch` |
| 129 | `BeginCombinatorialRelease` |
| 130 | `AdvanceCombinatorialRelease` |
| 131 | `InspectCombinatorialRelease` |
| 132 | `FinishCombinatorialRelease` |
| 133 | `ExportCombinatorialOutputs` |
| 134 | `SearchGoldenGateFusions` |
| 135 | `MaterializeGoldenGateFusions` |
| 136 | `ScoreLigationFidelity` |
| 137 | `InspectDomestication` |
| 138 | `ProposeDomestication` |
| 139 | `CheckDomestication` |
| 140 | `PrepareDomestication` |
| 141 | `ListLigationProfiles` |
| 142 | `ScorePublishedLigation` |
| 143 | `RankGoldenGateFusions` |
| 144 | `ScoreLigationCompetition` |
| 145 | `LoadRunEvidence` |
| 146 | `ProcessTopoActions` |
| 147 | `SimulateMultiSiteRecombination` |
| 148 | `ProjectMultiSiteRecombination` |
| 149 | `CheckMultiSiteRecombination` |
| 150 | `InspectCombinatorialMultiSite` |
| 151 | `PreviewTailingRecipe` |
| 152 | `InspectPcrDigest` |
| 153 | `InspectRunMaterial` |
| 154 | `InspectMaterialPcrDigest` |
| 155 | `InspectRunCatalogRestriction` |
| 156 | `BuildSubstitutedAssemblyIntent` |
| 157 | `InspectMismatchPcrDigest` |
| 158 | `BuildContributionEditIntent` |
| 159 | `DesignOverlapExtensionPrimers` |
| 160 | `InspectOverlapExtensionDigest` |
| 161 | `DesignEditedOverlapPrimers` |
| 162 | `PreviewHomologyRecipe` |
| 163 | `PreviewAnnealingRecipe` |
| 164 | `InspectTopoRun` |
| 165 | `InspectMaterialOrigin` |
| 166 | `MaterialPcrDocuments` |
| 167 | `PrepareMaterialPcrOutputs` |
| 168 | `ReviewSequenceExport` |
| 169 | `PrepareSequenceExport` |
| 170 | `ExportScfTrace` |
| 171 | `InstallMaterialDerivation` |
| 172 | `InspectMaterialDigest` |
| 173 | `PrepareSavedRunReview` |
| 174 | `RenderMapSvg` |
| 175 | `InspectMaterialMismatchPcr` |
| 176 | `InspectMaterialRepair` |
| 177 | `PreviewProtocolLayout` |
| 178 | `CaptureProtocolLayout` |
| 179 | `ExportProtocolLayout` |
| 180 | `PreviewProtocolQuantities` |
| 181 | `CaptureProtocolQuantities` |
| 182 | `OpenProtocolQuantities` |
| 183 | `ExportProtocolQuantities` |
| 184 | `ListProtocolQuantities` |
| 185 | `ListTopoGraphs` |
| 186 | `CaptureTopoGraph` |
| 187 | `InspectTopoGraph` |
| 188 | `InspectMaterialOverlapExtension` |
| 189 | `InspectMaterialTerminalProcessing` |
| 190 | `InspectMaterialAnnealedOligos` |

## Current result tags

These are indexes in the current Postcard enum. The current schema must change
when an internal enum layout changes.

| Tag | Result |
| ---: | --- |
| 0 | `Search` |
| 1 | `Restriction` |
| 2 | `Error` |
| 3 | `Pcr` |
| 4 | `Verification` |
| 5 | `DuplexDigest` |
| 6 | `Translation` |
| 7 | `CdsTranslation` |
| 8 | `Orfs` |
| 9 | `Alignment` |
| 10 | `PcrDocuments` |
| 11 | `PrimerEvidence` |
| 12 | `PrimerDesign` |
| 13 | `ProtocolPdf` |
| 14 | `GoldenGate` |
| 15 | `MismatchPcr` |
| 16 | `InversePcr` |
| 17 | `Gel` |
| 18 | `DiagnosticDigest` |
| 19 | `MismatchPcrPreview` |
| 20 | `WorkspaceSaved` |
| 21 | `WorkspaceLoaded` |
| 22 | `TraceValidation` |
| 23 | `FeatureLibraryScan` |
| 24 | `AnnotatedDocument` |
| 25 | `AnnotationImport` |
| 26 | `AnnotationExport` |
| 27 | `RepairPreview` |
| 28 | `AnnealPreview` |
| 29 | `PhysicalPreviewDocument` |
| 30 | `ViewTracks` |
| 31 | `DocumentSearch` |
| 32 | `TracksPrepared` |
| 33 | `AlignmentDocument` |
| 34 | `AlignmentExport` |
| 35 | `AlignmentRows` |
| 36 | `AlignmentSaved` |
| 37 | `ProteinProperties` |
| 38 | `TailPreview` |
| 39 | `OverlapProfile` |
| 40 | `RunCommit` |
| 41 | `RunReplay` |
| 42 | `CodonDesign` |
| 43 | `CodonCheck` |
| 44 | `DotPlot` |
| 45 | `Contigs` |
| 46 | `Methylation` |
| 47 | `CrisprGuides` |
| 48 | `CrisprOffTargets` |
| 49 | `EditTemplate` |
| 50 | `EditTemplateCheck` |
| 51 | `GenomeImportStarted` |
| 52 | `GenomeImportProgress` |
| 53 | `GenomeFastaFinished` |
| 54 | `GenomeRegistered` |
| 55 | `GenomeIntervalGuides` |
| 56 | `GenomeSearchState` |
| 57 | `GenomeSearchPage` |
| 58 | `GenomeDisposed` |
| 59 | `GenomeIntervalRead` |
| 60 | `RestrictionInspection` |
| 61 | `DigestProduct` |
| 62 | `DigestInstalled` |
| 63 | `CatalogRestrictionInspection` |
| 64 | `CatalogDigestProduct` |
| 65 | `CatalogDigestInstalled` |
| 66 | `RecombinationInspection` |
| 67 | `RecombinationTargetCheck` |
| 68 | `RecombinationProduct` |
| 69 | `ReadAlignmentsPrepared` |
| 70 | `ReadAlignmentPage` |
| 71 | `ReadAlignmentChecked` |
| 72 | `ReadAlignmentsExported` |
| 73 | `ReadAlignmentsDisposed` |
| 74 | `TopoSimulation` |
| 75 | `TopoOrientations` |
| 76 | `TopoProjection` |
| 77 | `TopoTarget` |
| 78 | `VariantsPrepared` |
| 79 | `VariantPage` |
| 80 | `VariantReferenceChecked` |
| 81 | `VariantAllelePreview` |
| 82 | `VariantsExported` |
| 83 | `VariantsDisposed` |
| 84 | `OverlapRecipePreview` |
| 85 | `LigationRecipePreview` |
| 86 | `AssemblyJunctionInspection` |
| 87 | `OverlapPrimersDesigned` |
| 88 | `AssemblyIntent` |
| 89 | `DigestFragmentSelected` |
| 90 | `PlanOligosExported` |
| 91 | `ProtocolHtml` |
| 92 | `CloningOutputsFinalized` |
| 93 | `GoldenGatePrimersDesigned` |
| 94 | `MutationIntent` |
| 95 | `ReadCheckCaptured` |
| 96 | `ReadCheckReplayed` |
| 97 | `ReadReviewPreview` |
| 98 | `ReadReviewCommitted` |
| 99 | `EditedAssemblyIntent` |
| 100 | `EditedAssemblyPrimers` |
| 101 | `SelectiveAssemblyPrimers` |
| 102 | `GoldenGateJunctionPrimers` |
| 103 | `PreparedGoldenGatePool` |
| 104 | `PreparedGoldenGateCandidate` |
| 105 | `GoldenGateReaction` |
| 106 | `GoldenGateReactionCandidate` |
| 107 | `PcrOutputs` |
| 108 | `MismatchPcrOutputs` |
| 109 | `ReviewedComparison` |
| 110 | `ReferenceFragment` |
| 111 | `ReferencePaste` |
| 112 | `TranslatedReference` |
| 113 | `TranslatedOutputs` |
| 114 | `PinnedGel` |
| 115 | `ExportedPairwiseComparison` |
| 116 | `SearchedPairwiseComparison` |
| 117 | `PreparedPairwiseContext` |
| 118 | `CombinatorialBatchStarted` |
| 119 | `CombinatorialBatchAdvanced` |
| 120 | `CombinatorialBatchCaptured` |
| 121 | `CombinatorialReleaseStarted` |
| 122 | `CombinatorialReleaseAdvanced` |
| 123 | `CombinatorialReleaseInspected` |
| 124 | `CombinatorialReleaseFinished` |
| 125 | `CombinatorialOutputsArchive` |
| 126 | `GoldenGateFusions` |
| 127 | `GoldenGateFusionMaterialized` |
| 128 | `LigationFidelity` |
| 129 | `DomesticationInspected` |
| 130 | `DomesticationProposed` |
| 131 | `DomesticationChecked` |
| 132 | `DomesticationPrepared` |
| 133 | `LigationProfiles` |
| 134 | `PublishedLigationScore` |
| 135 | `RankedGoldenGateFusions` |
| 136 | `LigationCompetition` |
| 137 | `HistoricalRunEvidence` |
| 138 | `TopoActionsProduct` |
| 139 | `MultiSiteRecombination` |
| 140 | `MultiSiteRecombinationProduct` |
| 141 | `MultiSiteRecombinationChecked` |
| 142 | `CombinatorialMultiSiteInspection` |
| 143 | `TailingRecipePreview` |
| 144 | `PcrDigestInspection` |
| 145 | `RunMaterialInspection` |
| 146 | `MaterialPcrDigestInspection` |
| 147 | `RunCatalogRestrictionInspection` |
| 148 | `SubstitutedAssemblyIntent` |
| 149 | `MismatchPcrDigestInspection` |
| 150 | `ContributionEditIntent` |
| 151 | `OverlapExtensionPrimersDesigned` |
| 152 | `OverlapExtensionDigestInspection` |
| 153 | `EditedOverlapPrimersDesigned` |
| 154 | `TopoRunInspection` |
| 155 | `MaterialOriginInspection` |
| 156 | `MaterialPcrDocuments` |
| 157 | `MaterialPcrOutputsPrepared` |
| 158 | `SequenceExportAssessment` |
| 159 | `SequenceExportPrepared` |
| 160 | `ScfTraceExported` |
| 161 | `MaterialDerivationInstalled` |
| 162 | `MaterialDigestInspection` |
| 163 | `SavedRunReviewPrepared` |
| 164 | `MapSvgExported` |
| 165 | `MaterialMismatchPcrPreview` |
| 166 | `MaterialRepairPreview` |
| 167 | `ProtocolLayoutPreview` |
| 168 | `ProtocolLayoutCaptured` |
| 169 | `ProtocolLayoutExported` |
| 170 | `ProtocolQuantitiesPreview` |
| 171 | `ProtocolQuantitiesCaptured` |
| 172 | `ProtocolQuantitiesOpened` |
| 173 | `ProtocolQuantitiesExported` |
| 174 | `ProtocolQuantitiesListed` |
| 175 | `TopoGraphsListed` |
| 176 | `TopoGraphCaptured` |
| 177 | `TopoGraphInspected` |
| 178 | `MaterialOverlapExtensionPreview` |
| 179 | `MaterialTerminalProcessingPreview` |
| 180 | `MaterialAnnealedOligosPreview` |

## Current source and evidence contracts

Exact PCR requests retain their explicit `template_strand`. Ordinary originals
use their supported reference frame; saved Run/Derived origins are freshly
reconstructed and checked against exact source events in the owning workspace.
PCR-product digest inspection preserves actual native termini and alternative
geometry. Serialized duplex summaries cannot establish template stock.

TOPO uses one chronological `Step` transcript for transfers, flap hydrolysis and
nick sealing. `SimulateTopo` and `InspectTopoRun::Simulate` inspect the full graph;
`ProcessTopoActions` uses the sole faithful tracked-material bridge.
Reference-only projection and target comparison remain distinct. Native graph
records retain complete source recipes and prefixes, including byproducts,
consumed activation and incomplete history. Graph-only continuation invents no
PCR input. See [graph transport](TOPO_GRAPH_TRANSPORT.md).

`LoadRunEvidence` returns the one `scientific-run-bundle-v1` attachment bundle.
Current physical graphs are retained by `physical_run_evidence`; producer stamps
and saved assessments remain documentary. Fresh reconstruction, not attachment
presence or historical status, establishes current scientific results.

Prepared ligation requests retain their explicit `CleavagePolicy`. Positional
checked Golden Gate preparation and candidate recutting use exact catalog enzyme
identity, actual modification states and bounded shared work. Unknown/impaired
source cleavage withholds the requested pool. Null original-generation
assumptions remain unknown; saved material cannot be redeclared as original.
Nested historical replay has its own bounded execution ledger plus the shared
ancestry limits; reported cleavage work describes the enclosing computation.

Plate assignments and reaction quantities are source-bound documentary records,
not scientific operations or caller-supplied readiness. Their preview/capture/open/
list/export receipts pin actual native sources and contexts. Captures preserve
biological heads and only add the exact reviewed attachment. Quantity arithmetic
completeness is separate from molecular verification; combined protocol exports
perform a fresh native check and retain both outcomes. No missing quantity or
reagent assumption is filled automatically.

Saved inverse-PCR proposals retain their original event and independent target.
An optional Run/Derived origin adds actual source reconstruction before PCR,
phosphorylation and intended ligation. Readonly response validation checks exact
correspondence; it does not replace that native execution.

## Related contracts

- [Retained run review](RUN_REVIEW.md)
- [Portable map export](MAP_SVG_EXPORT.md)
- [Ligation alternatives](LIGATION_ALTERNATIVES_EXPORT.md)
- [Native desktop IPC](DESKTOP_PROTOCOL.md)
- [Plate assignments](PROTOCOL_LAYOUT_TRANSPORT.md)
- [Combined quantity protocols](PROTOCOL_QUANTITIES_COMPOSITION.md)
- [Saved inverse-PCR sources](INVERSE_PCR_ORIGIN_TRANSPORT.md)

Compiler and scientific-fixture evidence is tracked in the implementation status.
Browser acceptance remains deferred until the requested code is complete.

`PrepareDomestication` carries an explicit nullable owning workspace (selector140
uses the workspace envelope). Its request carries a required-nullable native
origin. Original sources use null/null; saved Run/Derived sources require their
exact retained event and owner before preparation. Response validation checks the
caller-computed context and documentary correspondence without restoring material.
Library method values include explicit keyed TA/GC additions and current batch
records use schema2 with required-nullable typed end-family evidence.

## Target-free reusable overlap-extension PCR

`InspectMaterialOverlapExtension` (request188, result178) carries the exact owning
workspace, source event, explicit original/saved origin and six-oligo preparation.
Admission verifies that complete binding before sequence search. Its context binds
the source event, workspace, origin and preparation; the independently selected
outer product revision cannot be supplied as material authority. Installation uses
the existing generic material-derivation lane and freshly executes the retained
recipe. The inspector and saved derivation share one native OE executor.
Both branches and outer PCR require unique supported geometry. Inherited
incomplete evidence survives outer PCR and archive restart. No intended assembly
target is generated or assessed. Only the current worker schema is accepted.

## Reusable terminal processing

`InspectMaterialTerminalProcessing` (request189/result179) receives the exact
owning workspace, source event, original/saved origin and explicit terminal
action. The context includes the source sequence hash, exact event, origin hash
and action. Native kernels perform 5′ phosphate treatment or one base at both
3′ termini. The report retains source/result geometry and positional-modification
summaries; it does not supply material state. The current generic derivation
install lane independently repeats processing and rejects requests for new
primers. No automatic end repair, nick sealing, complement choice or independent
target is inserted.

## Role-bound reusable material

Current native material records and capsules use schema2, with explicit source
roles and complete role-bound source snapshots. Template operations require one
`template` source; supplied-oligo annealing requires `top_oligo`, then
`bottom_oligo`. The roles can refer to one self-complementary event. History
parents deduplicate the exact event identities; source roles retain both uses.
The shared ancestry/resource ledger and sharing closure visit every source.
There is no previous-format reader or singular-source fallback.

`InspectMaterialAnnealedOligos` (request190/result180) binds the owning workspace,
exact top/bottom events and explicit standard-oligo termini/search options. Its
selection retains a stable native candidate ID, preview context and actual
product revision. Preview and installation reuse native annealing and annotation
projection; installation independently repeats science and creates no primer
orders or invented assembly target. Alternatives and incomplete searches remain
conditional after selection and archive restart.
