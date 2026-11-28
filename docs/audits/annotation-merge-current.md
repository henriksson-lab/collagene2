# Non-destructive annotation merge: current owner audit

Source review, 28 November 2026. No Rust/fixture changes, builds, tests, browser or application requests. F19 currently has substantial append-and-review functionality; the missing capability is a unified, explicit merge decision with accumulated feature-level attribution, not a replacement annotation engine.

## Existing functionality

| Owner | Actual current behavior |
| --- | --- |
| `sequence/annotations/library.rs` | Versioned feature libraries retain exact source document revision, source feature ordinal/location/review state, entry and library identities. Scans distinguish concrete identity, compatible ambiguity, mismatches, complete whole-template coverage, repeats, palindromic orientation and work/hit incompleteness. Applying selected hits freshly reconstructs the complete scan, preserves existing annotations and DNA, marks additions for review and commits atomically. |
| Library feature qualifiers | New features already carry entry/library ID and revision, library name/version, source revision/ordinal, proposal ID and match counts/model. Attribution is per feature, not merely document metadata. Imported source identity remains a claim, not authenticated biological function. Document `feature_library:*` metadata records the last application; it is not a complete accumulated provenance ledger. |
| `formats/annotations/{mod,gtf}.rs` | GFF3/GTF/BED import has explicit exact-reference/seqid binding and retains source bytes, source hash, row ordinals, directives and source-record feature evidence. Selected annotations append without removing old features. Repeated exact Feature values are refused. Original source/application evidence accumulates by source-file hash in molecule metadata. Export validates retained row evidence against current geometry/qualifiers and refuses contradictory source evidence. |
| `documents/reference_edit.rs` | ApplyFeatureLibrary and ImportAnnotations are actual typed, replayable reference edits. Their native owners emit `AnnotationIdentity` during the same append execution: earlier feature ordinals remain retained and new ordinals are additions. F36 causal receipts already distinguish native annotation identity from merely similar labels/coordinates. |
| `sequence/annotations/batch.rs` and native editing | Explicit selected-feature operations preserve unselected annotations and DNA, with source-revision guards. Meaning/strand reassignment is distinct from retention; native receipts must continue to express that distinction. |
| Web and MCP | `web/viewer/sequence/feature_tools.rs` already offers library creation/import/export, bounded scan proposals and selected application, plus imported annotation preview/application. Worker-backed web edits retain native history. MCP sequence tools expose immutable library/dataset operations and add annotations into new documents, leaving originals intact. |

Existing duplicate definitions differ. Library application keys by feature kind, label and normalized resolved spans/orientations; it refuses an already-present key, including selected collisions within the same batch. Dataset application uses full Feature equality. Therefore two independently sourced but geometrically equivalent annotations can be refused by the library path or appended separately by the file path. Neither behavior silently overwrites existing features, but neither offers an explicit evidence-collation decision.

## Actual missing capability

1. **Explicit, reviewable conflict/duplicate decisions.** There is no common preview showing exact duplicate, same geometry with different meanings/qualifiers, independent attribution, orientation disagreement and merely overlapping features as different cases. Overlap alone must not collapse two features. Identical labels must not establish feature identity.
2. **Idempotent non-destructive application.** Selecting one already-present library proposal currently rejects the whole selection, including unrelated new features. A declared Skip decision can provide a useful atomic merge without claiming new evidence or silently erasing annotations. It must be recorded explicitly, not reinterpret an error as success.
3. **Accumulated provenance on the resulting feature.** Existing library and imported-file attribution is real, but there is no common bounded owner for multiple accepted source contributions to one annotation. Replacing the document's last library application metadata does not preserve an accumulated feature-level history by itself.
4. **Annotation-only transfer between existing documents.** Current library extraction/re-scan or file export/import can transfer annotations indirectly. There is no direct selected-feature merge from an independently bound existing document, with a declared exact coordinate correspondence and no DNA replacement. PasteReference is an actual DNA edit and must not become a disguised annotation merge.
5. **Unified UI/MCP review and provenance inspection.** Existing controls handle their respective proposal families. They do not present one explicit conflict policy and inspect multiple retained contributions to a selected feature.

## Proposed cohesive extension

Extend `sequence::annotations` with a `merge` owner, reusing current library scan/import proposal validation and model projection/editing. Keep one typed ReferenceEdit operation for the final merge and one current merge request/report representation across native/worker/web/MCP. Do not synthesize a DesignPlan, a second history or another sequence matcher.

A concrete initial contract can be:

- `Request { target_revision, source: SourceEvidence, selections: Vec<Selection>, limits }`. SourceEvidence identifies one actual library+scan, imported dataset, or explicitly bound donor document. Selection pins an actual proposal ID/source feature ordinal and a decision: Append, Skip, or AddProvenance to an exact existing target feature ordinal plus its content revision. Mixed-source batches can follow only with a single shared native work/retention ledger.
- `preview(target, source, request) -> MergeReview` freshly validates the current source owner. Each candidate retains source identity, exact proposed feature, matched target ordinals/content pins and a typed relationship. Distinguish exact whole-feature duplicate, exact geometry with qualifier/type disagreement, biological-path uncertainty and overlap. Preview is evidence, not permission to overwrite or an authenticated function statement.
- `apply(target, source, request) -> AppliedMerge { molecule, receipt }` recomputes the same source evidence and executes the explicit decisions atomically. Existing feature order/content and every DNA call remain unchanged for Append/Skip. AddProvenance changes attribution only and must not overwrite biological qualifiers or clear needs_review. A semantic merge of qualifier values should be a separately explicit later decision with actual review-required annotation reassignment.
- Receipt records actual original/destination ordinals and each append/skip/attribution action, emitted in the same execution. Reuse `AnnotationIdentity`; preserve previous identity for genuinely attribution-only additions, while explicit semantic redefinitions sever it. Do not infer retained identity after comparing labels, hashes or final arrays.

For direct donor transfer, begin with an explicit full-record identity/reverse-complement/rotation correspondence checked against actual DNA and topology, then use canonical location projection. The requested transform is part of the declared map; repeated equivalent origins are not automatically unique biological correspondence. Partial homology, indel alignment and uncertain remapping remain unsupported until genuinely covered by an existing canonical owner. Joined/reverse/origin-spanning features require actual biological paths and projection receipts; fuzzy/unresolved locations cannot be reduced to a broad min/max interval.

## Provenance storage and bounds

Choose the canonical provenance owner before implementation. A typed model feature field would require one direct current molecule/schema migration and updating all constructors/codecs. A reserved feature-qualifier envelope can avoid a model dependency cycle, but must have one typed parser/validator and one current encoding, not parallel ad-hoc strings and native records. Extend/update current internal library attribution directly if superseded; existing external GFF/GTF/BED source-row evidence remains external-format functionality, not a legacy internal decoder. Do not copy entire library/dataset/file bytes into every feature: retain bounded content-addressed source references and exact contributing ordinals; shared original evidence remains in the native workspace/document source owner.

Before cloning, admit target/source model trees, selected count, comparison candidates, aggregate locations, attribution count/text/bytes and final output size. Reuse library work200M/hits10k and file8MiB/rows10k guards, but use one aggregate merge allowance rather than multiplying full fresh scans per selection. A single source scan/import must be reconstructed once per application, then all decisions derived from that evidence. Incomplete scans cannot certify absence or prove that an apparent duplicate is the only corresponding feature. Provenance is documentary attribution; imported qualifiers never promote review status or confer native event ancestry.

## Independent fixtures needed

- Existing features and DNA survive a mixed Append+explicit Skip merge exactly; a late invalid decision commits nothing.
- Two features with the same label at different positions remain separate; overlapping promoters/CDSs remain separate; same geometry with conflicting qualifiers requires an explicit decision.
- Two independently bound identical-DNA donors remain separately attributed; source/document/feature ordinal changes invalidate captured decisions.
- Reapplying an exact source contribution is idempotent under explicit Skip; a different independent contribution can accumulate without deleting either source record or copying source bytes per feature.
- Reverse, joined and origin-crossing transfers use literal expected coordinates; unresolved/fuzzy correspondence is explicit uncertainty/refusal.
- Ambiguous library matches, incomplete scan budgets and palindromic orientations retain their current uncertainty; no automatic review clearance or functional validation.
- Native merge receipts retain existing ordinals and record actual additions; later semantic reassignment does not masquerade as inherited annotation identity in F36.
- Current workspace/candidate/history serialization, immutable MCP application and web undo/reopen retain accumulated attribution. External exports either retain supported provenance or report declared loss.

The initial product should deliver conflict-aware append/skip and attribution inspection before attempting generic qualifier reconciliation. Most scanning, import/export, source binding and append infrastructure already exists and should be extended directly.
