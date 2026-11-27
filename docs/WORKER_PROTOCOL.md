# Binary worker protocol

The web worker and native adapters share Rust DTOs in `core::protocol` and the
same deterministic scientific engines. Current schemas are worker **96**,
candidate **25** and native workspace plan **23**; the verifier is **science-26**.
Only current internal encodings are accepted. Superseded development formats have
no fallback decoder, migration or capability ladder.

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

These are indexes in the current Postcard enum, not a compatibility promise for
superseded development schemas. Change the current schema when changing layouts.

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
| 86 | `ProjectAnnotatedTopo` |
| 87 | `CheckAnnotatedTopoTarget` |
| 88 | `PrepareVariants` |
| 89 | `VariantPage` |
| 90 | `CheckVariantReference` |
| 91 | `PreviewVariantAllele` |
| 92 | `ExportVariants` |
| 93 | `DisposeVariants` |
| 94 | `PreviewOverlapRecipe` |
| 95 | `PreviewLigationRecipe` |
| 96 | `InspectAssemblyJunctions` |
| 97 | `DesignOverlapPrimers` |
| 98 | `BuildAssemblyIntent` |
| 99 | `SelectDigestFragment` |
| 100 | `ExportPlanOligos` |
| 101 | `ExportProtocolHtml` |
| 102 | `FinalizeCloningOutputs` |
| 103 | `DesignGoldenGatePrimers` |
| 104 | `BuildMutationIntent` |
| 105 | `CaptureReadCheck` |
| 106 | `ReplayReadCheck` |
| 107 | `PreviewReadReview` |
| 108 | `CommitReadReview` |
| 109 | `BuildEditedAssemblyIntent` |
| 110 | `DesignEditedAssemblyPrimers` |
| 111 | `DesignSelectiveAssemblyPrimers` |
| 112 | `DesignGoldenGateJunctionPrimers` |
| 113 | `InspectPreparedGoldenGatePool` |
| 114 | `SelectPreparedGoldenGateCandidate` |
| 115 | `InspectGoldenGateReaction` |
| 116 | `SelectGoldenGateReaction` |
| 117 | `PreparePcrOutputs` |
| 118 | `PrepareMismatchPcrOutputs` |
| 119 | `CompareSequences` |
| 120 | `CopyReferenceFragment` |
| 121 | `PasteReferenceFragment` |
| 122 | `TranslateReference` |
| 123 | `PrepareTranslatedReference` |
| 124 | `PredictPinnedGel` |
| 125 | `ExportPairwiseComparison` |
| 126 | `SearchPairwiseComparison` |
| 127 | `PreparePairwiseContext` |
| 128 | `BeginCombinatorialBatch` |
| 129 | `AdvanceCombinatorialBatch` |
| 130 | `CaptureCombinatorialBatch` |
| 131 | `BeginCombinatorialRelease` |
| 132 | `AdvanceCombinatorialRelease` |
| 133 | `InspectCombinatorialRelease` |
| 134 | `FinishCombinatorialRelease` |
| 135 | `ExportCombinatorialOutputs` |
| 136 | `SearchGoldenGateFusions` |
| 137 | `MaterializeGoldenGateFusions` |
| 138 | `ScoreLigationFidelity` |
| 139 | `InspectDomestication` |
| 140 | `ProposeDomestication` |
| 141 | `CheckDomestication` |
| 142 | `PrepareDomestication` |
| 143 | `ListLigationProfiles` |
| 144 | `ScorePublishedLigation` |
| 145 | `RankGoldenGateFusions` |
| 146 | `ScoreLigationCompetition` |
| 147 | `SealTopoNicks` |
| 148 | `ProjectProcessedTopo` |
| 149 | `LoadProcessedRunEvidence` |
| 150 | `ProcessTopoActions` |
| 151 | `InspectTopoActions` |
| 152 | `SimulateMultiSiteRecombination` |
| 153 | `ProjectMultiSiteRecombination` |
| 154 | `CheckMultiSiteRecombination` |
| 155 | `InspectCombinatorialMultiSite` |
| 156 | `PreviewTailingRecipe` |
| 157 | `InspectPcrDigest` |
| 158 | `InspectRunMaterial` |
| 159 | `InspectRunPcrDigest` |
| 160 | `InspectRunCatalogRestriction` |
| 161 | `BuildSubstitutedAssemblyIntent` |
| 162 | `InspectMismatchPcrDigest` |
| 163 | `BuildContributionEditIntent` |
| 164 | `DesignOverlapExtensionPrimers` |
| 165 | `InspectOverlapExtensionDigest` |
| 166 | `DesignEditedOverlapPrimers` |
| 167 | `PreviewHomologyRecipe` |
| 168 | `PreviewAnnealingRecipe` |
| 169 | `InspectTopoRun` |
| 170 | `InspectMaterialOrigin` |
| 171 | `MaterialPcrDocuments` |
| 172 | `PrepareMaterialPcrOutputs` |
| 173 | `ReviewSequenceExport` |
| 174 | `PrepareSequenceExport` |
| 175 | `ExportScfTrace` |
| 176 | `InstallMaterialDerivation` |
| 177 | `InspectMaterialDigest` |
| 178 | `PrepareSavedRunReview` |
| 179 | `RenderMapSvg` |
| 180 | `InspectMaterialMismatchPcr` |
| 181 | `InspectMaterialRepair` |

## Current result tags

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
| 78 | `AnnotatedTopoProjection` |
| 79 | `AnnotatedTopoTarget` |
| 80 | `VariantsPrepared` |
| 81 | `VariantPage` |
| 82 | `VariantReferenceChecked` |
| 83 | `VariantAllelePreview` |
| 84 | `VariantsExported` |
| 85 | `VariantsDisposed` |
| 86 | `OverlapRecipePreview` |
| 87 | `LigationRecipePreview` |
| 88 | `AssemblyJunctionInspection` |
| 89 | `OverlapPrimersDesigned` |
| 90 | `AssemblyIntent` |
| 91 | `DigestFragmentSelected` |
| 92 | `PlanOligosExported` |
| 93 | `ProtocolHtml` |
| 94 | `CloningOutputsFinalized` |
| 95 | `GoldenGatePrimersDesigned` |
| 96 | `MutationIntent` |
| 97 | `ReadCheckCaptured` |
| 98 | `ReadCheckReplayed` |
| 99 | `ReadReviewPreview` |
| 100 | `ReadReviewCommitted` |
| 101 | `EditedAssemblyIntent` |
| 102 | `EditedAssemblyPrimers` |
| 103 | `SelectiveAssemblyPrimers` |
| 104 | `GoldenGateJunctionPrimers` |
| 105 | `PreparedGoldenGatePool` |
| 106 | `PreparedGoldenGateCandidate` |
| 107 | `GoldenGateReaction` |
| 108 | `GoldenGateReactionCandidate` |
| 109 | `PcrOutputs` |
| 110 | `MismatchPcrOutputs` |
| 111 | `ReviewedComparison` |
| 112 | `ReferenceFragment` |
| 113 | `ReferencePaste` |
| 114 | `TranslatedReference` |
| 115 | `TranslatedOutputs` |
| 116 | `PinnedGel` |
| 117 | `ExportedPairwiseComparison` |
| 118 | `SearchedPairwiseComparison` |
| 119 | `PreparedPairwiseContext` |
| 120 | `CombinatorialBatchStarted` |
| 121 | `CombinatorialBatchAdvanced` |
| 122 | `CombinatorialBatchCaptured` |
| 123 | `CombinatorialReleaseStarted` |
| 124 | `CombinatorialReleaseAdvanced` |
| 125 | `CombinatorialReleaseInspected` |
| 126 | `CombinatorialReleaseFinished` |
| 127 | `CombinatorialOutputsArchive` |
| 128 | `GoldenGateFusions` |
| 129 | `GoldenGateFusionMaterialized` |
| 130 | `LigationFidelity` |
| 131 | `DomesticationInspected` |
| 132 | `DomesticationProposed` |
| 133 | `DomesticationChecked` |
| 134 | `DomesticationPrepared` |
| 135 | `LigationProfiles` |
| 136 | `PublishedLigationScore` |
| 137 | `RankedGoldenGateFusions` |
| 138 | `LigationCompetition` |
| 139 | `TopoNicksSealed` |
| 140 | `ProcessedTopoProduct` |
| 141 | `HistoricalProcessedRunEvidence` |
| 142 | `TopoActionsProduct` |
| 143 | `TopoActionsInspected` |
| 144 | `MultiSiteRecombination` |
| 145 | `MultiSiteRecombinationProduct` |
| 146 | `MultiSiteRecombinationChecked` |
| 147 | `CombinatorialMultiSiteInspection` |
| 148 | `TailingRecipePreview` |
| 149 | `PcrDigestInspection` |
| 150 | `RunMaterialInspection` |
| 151 | `RunPcrDigestInspection` |
| 152 | `RunCatalogRestrictionInspection` |
| 153 | `SubstitutedAssemblyIntent` |
| 154 | `MismatchPcrDigestInspection` |
| 155 | `ContributionEditIntent` |
| 156 | `OverlapExtensionPrimersDesigned` |
| 157 | `OverlapExtensionDigestInspection` |
| 158 | `EditedOverlapPrimersDesigned` |
| 159 | `TopoRunInspection` |
| 160 | `MaterialOriginInspection` |
| 161 | `MaterialPcrDocuments` |
| 162 | `MaterialPcrOutputsPrepared` |
| 163 | `SequenceExportAssessment` |
| 164 | `SequenceExportPrepared` |
| 165 | `ScfTraceExported` |
| 166 | `MaterialDerivationInstalled` |
| 167 | `MaterialDigestInspection` |
| 168 | `SavedRunReviewPrepared` |
| 169 | `MapSvgExported` |
| 170 | `MaterialMismatchPcrPreview` |
| 171 | `MaterialRepairPreview` |

## Related contracts

- [Retained run review](RUN_REVIEW.md)
- [Portable map export](MAP_SVG_EXPORT.md)
- [Ligation alternatives](LIGATION_ALTERNATIVES_EXPORT.md)
- [Native desktop IPC](DESKTOP_PROTOCOL.md)

Compiler and scientific-fixture evidence is tracked in the implementation status.
Browser acceptance remains deferred until the requested code is complete.

Mismatch-PCR and repair previews (jobs180/181, results170/171) use the explicit
128 MiB compound lane because their source origins may close native ancestry.
They bind exact source/origin and operation settings; native replay establishes
material at installation. Map SVG (job179/result169) retains the ordinary limit.
