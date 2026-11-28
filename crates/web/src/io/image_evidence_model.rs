//! Documentary-image capture and workspace CAS pins. No image codec runs here.
use collagene_core::{
    documents::{
        image_evidence::{Binding, MAX_SOURCE_BYTES, Prepared, Summary},
        workspace::Workspace,
    },
    material::lineage::workspace_fingerprint,
};
use std::sync::{
    Arc,
    atomic::{AtomicBool, Ordering},
};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Context {
    pub epoch: u64,
    pub workspace_fingerprint: String,
    pub association: Binding,
    pub document_name: String,
}
impl Context {
    pub fn new(workspace: &Workspace, epoch: u64) -> Result<Self, String> {
        let fingerprint = workspace_fingerprint(workspace).map_err(|e| e.to_string())?;
        let document = workspace
            .documents
            .iter()
            .find(|d| d.id == workspace.active_document)
            .ok_or("The active document is unavailable.")?;
        let source = workspace.head(&document.id).map_err(|e| e.to_string())?;
        Ok(Self {
            epoch,
            workspace_fingerprint: fingerprint,
            association: Binding {
                document_id: document.id.clone(),
                revision: document.head.clone(),
            },
            document_name: source.name.clone(),
        })
    }
    pub fn matches(&self, current: &Self) -> Result<(), String> {
        if self != current {
            return Err("The workspace, document or retained attachments changed. Close and reopen image evidence before continuing.".into());
        }
        Ok(())
    }
}
pub struct Snapshot {
    pub workspace: Workspace,
    pub context: Context,
}
impl Snapshot {
    pub fn new(workspace: Workspace, epoch: u64) -> Result<Self, String> {
        let context = Context::new(&workspace, epoch)?;
        Ok(Self { workspace, context })
    }
}
#[derive(Clone, Debug)]
pub struct FileCapture {
    pub context: Context,
    pub name: String,
    pub original_hash: String,
    pub source_bytes: usize,
    pub revision: String,
}
impl FileCapture {
    pub fn new(
        context: Context,
        name: String,
        original_hash: String,
        source_bytes: usize,
    ) -> Result<Self, String> {
        if source_bytes == 0 || source_bytes > MAX_SOURCE_BYTES {
            return Err("Choose a nonempty image file of at most 16 MiB.".into());
        }
        let revision = context
            .association
            .preparation_revision(&name, &original_hash)
            .map_err(|e| e.to_string())?;
        Ok(Self {
            context,
            name,
            original_hash,
            source_bytes,
            revision,
        })
    }
    /// A private worker's result must identify the actual file and selected document,
    /// not just equal pixels or a matching filename. Prepared::summary checks asset hashes.
    pub fn admit(&self, prepared: &Prepared) -> Result<Summary, String> {
        let summary = prepared.summary().map_err(|e| e.to_string())?;
        if summary.record.name != self.name
            || summary.record.association != self.context.association
            || summary.record.original.hash != self.original_hash
            || summary.record.original.bytes != self.source_bytes as u64
            || summary
                .record
                .association
                .preparation_revision(&summary.record.name, &summary.record.original.hash)
                .map_err(|e| e.to_string())?
                != self.revision
        {
            return Err("The prepared image belongs to a different file, name or document.".into());
        }
        Ok(summary)
    }
}
pub struct Install {
    pub context: Context,
    pub prepared: Arc<Prepared>,
    pub alive: Arc<AtomicBool>,
}
impl Install {
    pub fn ensure_live(&self) -> Result<(), String> {
        if !self.alive.load(Ordering::SeqCst) {
            return Err("The image evidence dialog closed before attachment admission.".into());
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::documents::workspace::WorkspaceAttachment;
    use collagene_model::{Molecule, Topology};
    fn frame() -> Workspace {
        Workspace::from_documents(
            "Image evidence capture",
            vec![
                Molecule::dna("first", "ACGTACGT", Topology::Circular).unwrap(),
                Molecule::dna("second", "AGTCAGTC", Topology::Linear).unwrap(),
            ],
            0,
        )
        .unwrap()
    }
    #[test]
    fn workspace_epoch_active_document_and_unrelated_asset_edits_invalidate_image_cas() {
        let original = frame();
        let context = Context::new(&original, 7).unwrap();
        context
            .matches(&Context::new(&original, 7).unwrap())
            .unwrap();
        assert!(
            context
                .matches(&Context::new(&original, 8).unwrap())
                .is_err()
        );
        let mut changed = original.clone();
        changed.active_document = changed.documents[1].id.clone();
        assert!(
            context
                .matches(&Context::new(&changed, 7).unwrap())
                .is_err()
        );
        let mut changed = original.clone();
        changed.name.push_str(" renamed");
        assert!(
            context
                .matches(&Context::new(&changed, 7).unwrap())
                .is_err()
        );
        let mut changed = original;
        let bytes = b"unrelated attachment edit".to_vec();
        let hash = blake3::hash(&bytes).to_hex().to_string();
        changed.attachments.insert(
            hash,
            WorkspaceAttachment {
                kind: "documentary-note".into(),
                name: "note.txt".into(),
                media_type: "text/plain".into(),
                bytes,
            },
        );
        assert!(
            context
                .matches(&Context::new(&changed, 7).unwrap())
                .is_err(),
            "a pending image must not overwrite another attachment edit"
        );
    }
    #[test]
    fn image_admission_binds_actual_file_bytes_name_size_and_selected_document() {
        const SOURCE: &[u8] = &[
            137, 80, 78, 71, 13, 10, 26, 10, 0, 0, 0, 13, 73, 72, 68, 82, 0, 0, 0, 2, 0, 0, 0, 1,
            8, 2, 0, 0, 0, 123, 64, 232, 221, 0, 0, 0, 15, 73, 68, 65, 84, 120, 156, 99, 248, 207,
            192, 192, 240, 159, 1, 0, 7, 255, 1, 255, 1, 127, 137, 167, 0, 0, 0, 0, 73, 69, 78, 68,
            174, 66, 96, 130,
        ];
        const ANOTHER_SOURCE: &[u8] = &[
            137, 80, 78, 71, 13, 10, 26, 10, 0, 0, 0, 13, 73, 72, 68, 82, 0, 0, 0, 2, 0, 0, 0, 1,
            8, 2, 0, 0, 0, 123, 64, 232, 221, 0, 0, 0, 37, 116, 69, 88, 116, 78, 111, 116, 101, 0,
            115, 97, 109, 101, 32, 112, 105, 120, 101, 108, 115, 44, 32, 97, 110, 111, 116, 104,
            101, 114, 32, 115, 111, 117, 114, 99, 101, 32, 102, 105, 108, 101, 90, 143, 117, 214,
            0, 0, 0, 15, 73, 68, 65, 84, 120, 156, 99, 248, 207, 192, 192, 240, 159, 1, 0, 7, 255,
            1, 255, 1, 127, 137, 167, 0, 0, 0, 0, 73, 69, 78, 68, 174, 66, 96, 130,
        ];
        let context = Context::new(&frame(), 9).unwrap();
        let hash = blake3::hash(SOURCE).to_hex().to_string();
        let capture = FileCapture::new(
            context.clone(),
            "capture.png".into(),
            hash.clone(),
            SOURCE.len(),
        )
        .unwrap();
        assert_eq!(
            capture.revision,
            context
                .association
                .preparation_revision("capture.png", &hash)
                .unwrap()
        );
        // Fixture construction imitates the worker; UI admission itself has no codec call.
        let prepared =
            Prepared::new("capture.png", context.association.clone(), SOURCE.to_vec()).unwrap();
        let summary = capture.admit(&prepared).unwrap();
        assert_eq!(
            (
                summary.record.original.width,
                summary.record.original.height
            ),
            (2, 1)
        );
        prepared
            .preview()
            .unwrap()
            .validate_binding(&summary)
            .unwrap();
        let another = Prepared::new(
            "capture.png",
            context.association.clone(),
            ANOTHER_SOURCE.to_vec(),
        )
        .unwrap();
        assert!(
            capture.admit(&another).is_err(),
            "equal displayed pixels are not the same source file"
        );
        let renamed = FileCapture::new(
            context.clone(),
            "another label".into(),
            hash.clone(),
            SOURCE.len(),
        )
        .unwrap();
        assert!(renamed.admit(&prepared).is_err());
        let mut another_context = context.clone();
        another_context.association.document_id = "another-document".into();
        let wrong_owner = FileCapture::new(
            another_context,
            "capture.png".into(),
            hash.clone(),
            SOURCE.len(),
        )
        .unwrap();
        assert!(wrong_owner.admit(&prepared).is_err());
        let wrong_size = FileCapture::new(
            context.clone(),
            "capture.png".into(),
            hash.clone(),
            SOURCE.len() + 1,
        )
        .unwrap();
        assert!(wrong_size.admit(&prepared).is_err());
        let alive = Arc::new(AtomicBool::new(true));
        let install = Install {
            context: context.clone(),
            prepared: Arc::new(prepared),
            alive: alive.clone(),
        };
        install.ensure_live().unwrap();
        alive.store(false, Ordering::SeqCst);
        assert!(install.ensure_live().is_err());
        assert!(FileCapture::new(context.clone(), "capture.png".into(), hash.clone(), 0).is_err());
        assert!(
            FileCapture::new(context, "capture.png".into(), hash, MAX_SOURCE_BYTES + 1).is_err()
        );
    }
}
