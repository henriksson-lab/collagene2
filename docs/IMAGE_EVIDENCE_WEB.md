# Documentary images

Open a sequence document and choose **File → Images…**. Choose a JPEG, PNG or TIFF file up to 16 MiB, give it a name, and prepare its preview. Review the preview before choosing **Attach image**. The image is associated with the selected document revision.

Images are retained in native workspace (`.cg2design`) files and local workspace autosave. Standalone sequence exports, including `.cg2`, GenBank and FASTA, do not contain these attachments. Before attachment, a private persistence worker verifies that the complete proposed workspace can be encoded within the native archive limits. The application then checks that the workspace, selected document and dialog are still current, installs the attachments, and uses the existing generation-guarded autosave lane.

Reopen **File → Images…** to see images associated with the active document, including retained historical revisions. **Preview** decodes the retained original in a private worker and verifies the resulting preview. **Download original** saves the exact retained source bytes, including metadata and any additional TIFF pages. Its normal browser or native save destination is used.

Previews are documentary, uncalibrated display images. They do not verify cloning results or compare observed and predicted gels. Previews may be downsampled or reduced in sample precision; only the first TIFF page is shown. The dialog reports applicable transformations, and the original remains unchanged.

## Implementation contract

- Source decoding and preview generation run in isolated workers. The interface hashes the bounded file bytes and validates source identity, metadata and preview binding; it does not decode JPEG/TIFF or decompress archives on the main thread.
- Preparation captures the file hash and size, display name, document ID and exact revision. Worker replies must match their request ID and canonical revision context before they can become previews or attachment candidates.
- A fresh prepared response supplies its checked preview without another codec pass. An archived record is reopened through the worker before display.
- Attachment captures the complete workspace fingerprint and workspace epoch. Archive encoding happens before any editor mutation. After the worker returns, the fingerprint, epoch, active document and dialog lifetime are checked again. Only attachments are installed; sequence data, feature origins, document histories and cloning candidates remain owned by their existing workflows.
- Closing a dialog invalidates pending codec results and prevents pending attachments from being admitted. Exact-original downloads use the existing guarded export path; an acknowledged native save may finish after the dialog closes.
- Preview object URLs are scoped to their displayed preview and revoked when it is released.

Native/WASM checks and302focused cases pass; the [integration record](audits/image-evidence-current-integration.json) lists exact scope and pending publication CI. Browser and installed-GUI review remain deferred.
