//! Pure source/preview pins for explicit reusable PCR/digest installation.
//! This module never edits a workspace, writes recovery storage or reports commit.
use collagene_core::{
    Molecule,
    analysis::Primer,
    catalog_digest::CatalogSelection,
    material_derivation::{self as native, Derivation, PcrSelection, SourceOrigin},
    restriction_inspection::{FragmentSummary, InspectionOptions},
    run_lineage::EventBinding,
    run_material::MaterialOrigin,
    verification::PcrDocument,
    workspace::Workspace,
};
use serde::Serialize;
use std::{io::Write, sync::Arc};

const MAX_PIN_BYTES: usize = 128 * 1024 * 1024;
fn hash(value: &impl Serialize) -> Result<String, String> {
    struct Bound {
        bytes: usize,
        hash: blake3::Hasher,
    }
    impl Write for Bound {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.bytes = self
                .bytes
                .checked_add(bytes.len())
                .filter(|n| *n <= MAX_PIN_BYTES)
                .ok_or_else(|| std::io::Error::other("derivation pin exceeds bounded context"))?;
            self.hash.update(bytes);
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut bound = Bound {
        bytes: 0,
        hash: blake3::Hasher::new(),
    };
    serde_json::to_writer(&mut bound, value).map_err(|e| e.to_string())?;
    Ok(bound.hash.finalize().to_hex().to_string())
}
fn hex(text: &str) -> bool {
    text.len() == 64
        && text
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn identifier(text: &str) -> bool {
    !text.is_empty() && text.len() <= 512 && !text.chars().any(char::is_control)
}
#[derive(Debug, Clone)]
pub struct SourcePin {
    event: EventBinding,
    epoch: u64,
    origin: Option<Arc<MaterialOrigin>>,
    origin_revision: Option<String>,
}
impl SourcePin {
    pub fn capture(
        event: EventBinding,
        epoch: u64,
        source: &Molecule,
        origin: Option<Arc<MaterialOrigin>>,
    ) -> Result<Self, String> {
        collagene_core::model::bounds::model_preflight(&[source]).map_err(|e| e.to_string())?;
        source.validate().map_err(|e| e.to_string())?;
        if !identifier(&event.document_id)
            || !identifier(&event.event_id)
            || !hex(&event.revision)
            || event.revision != source.revision()
        {
            return Err("Choose the exact source document/event/revision before opening a reusable product.".into());
        }
        if let Some(origin) = &origin {
            origin.validate_documentary().map_err(|e| e.to_string())?;
            if origin.event() != &event {
                return Err("Saved material belongs to another exact source event.".into());
            }
        }
        let origin_revision = origin.as_ref().map(|o| hash(o.as_ref())).transpose()?;
        Ok(Self {
            event,
            epoch,
            origin,
            origin_revision,
        })
    }
    pub fn event(&self) -> &EventBinding {
        &self.event
    }
    pub fn epoch(&self) -> u64 {
        self.epoch
    }
    pub fn origin(&self) -> Option<&MaterialOrigin> {
        self.origin.as_deref()
    }
    pub fn validate_current(
        &self,
        event: &EventBinding,
        epoch: u64,
        source: &Molecule,
        origin: Option<&MaterialOrigin>,
    ) -> Result<(), String> {
        if event != &self.event
            || epoch != self.epoch
            || source.revision() != self.event.revision
            || origin.map(hash).transpose()? != self.origin_revision
        {
            return Err("The source event, workspace, or saved-material origin changed; recompute and select the product again.".into());
        }
        Ok(())
    }
    pub(crate) fn native_origin(&self) -> SourceOrigin {
        match self.origin.as_deref() {
            Some(origin) => SourceOrigin::Saved(origin.clone()),
            None => SourceOrigin::OriginalReference,
        }
    }
    pub(crate) fn attest(&self, workspace: &Workspace) -> Result<(), String> {
        if let Some(origin) = self.origin.as_deref() {
            origin
                .validate_for_workspace(workspace, &self.event)
                .map_err(|e| e.to_string())?;
        }
        let source = workspace
            .revisions
            .get(&self.event.revision)
            .ok_or("The captured source revision is no longer retained.")?;
        self.validate_current(&self.event, self.epoch, source, self.origin.as_deref())
    }
}
#[derive(Debug, Clone)]
pub struct PcrChoice {
    source: SourcePin,
    forward: Primer,
    reverse: Primer,
    selected_index: usize,
    preview: Arc<PcrDocument>,
    preview_revision: String,
    preview_context: String,
}
impl PcrChoice {
    pub fn capture(
        source: SourcePin,
        forward: Primer,
        reverse: Primer,
        selected_index: usize,
        preview: Arc<PcrDocument>,
    ) -> Result<Self, String> {
        if selected_index >= 10_000 {
            return Err("The selected PCR geometry exceeds the bounded preview profile.".into());
        }
        preview.molecule.validate().map_err(|e| e.to_string())?;
        let preview_revision = preview.molecule.revision();
        let preview_context = hash(&(&forward, &reverse, selected_index, preview.as_ref()))?;
        Ok(Self {
            source,
            forward,
            reverse,
            selected_index,
            preview,
            preview_revision,
            preview_context,
        })
    }
    pub fn source(&self) -> &SourcePin {
        &self.source
    }
    pub fn preview_revision(&self) -> &str {
        &self.preview_revision
    }
    pub fn validate_current(
        &self,
        event: &EventBinding,
        epoch: u64,
        source: &Molecule,
        origin: Option<&MaterialOrigin>,
        forward: &Primer,
        reverse: &Primer,
        selected_index: usize,
        preview: &PcrDocument,
    ) -> Result<(), String> {
        self.source.validate_current(event, epoch, source, origin)?;
        if hash(&(forward, reverse, selected_index, preview))? != self.preview_context {
            return Err("PCR oligos, selected geometry, or preview changed; recompute before opening reusable material.".into());
        }
        Ok(())
    }
    pub fn request(
        &self,
        workspace: &Workspace,
        name: &str,
    ) -> Result<PendingInstallation, String> {
        self.request_with_primers(workspace, name, false)
    }
    pub fn request_with_primers(
        &self,
        workspace: &Workspace,
        name: &str,
        save_primers: bool,
    ) -> Result<PendingInstallation, String> {
        self.source.attest(workspace)?;
        PendingInstallation::new(
            workspace,
            &self.source,
            native::Request {
                source: self.source.event.clone(),
                origin: self.source.native_origin(),
                name: name.into(),
                save_primers,
                derivation: Derivation::ExactPcr {
                    forward: self.forward.clone(),
                    reverse: self.reverse.clone(),
                    selection: PcrSelection {
                        selected_index: self.selected_index,
                        expected_preview_revision: self.preview_revision.clone(),
                        geometry: self.preview.geometry.clone(),
                    },
                },
            },
        )
    }
}
#[derive(Debug, Clone)]
pub struct DigestChoice {
    source: SourcePin,
    selection: CatalogSelection,
    options: InspectionOptions,
    fragment: FragmentSummary,
    inspection_revision: String,
    preview_context: String,
}
impl DigestChoice {
    /// The caller supplies a validated, current native inspection context. This
    /// records its choice only; native installation recomputes actual cleavage.
    pub fn capture(
        source: SourcePin,
        selection: CatalogSelection,
        options: InspectionOptions,
        inspection_revision: String,
        fragment: FragmentSummary,
    ) -> Result<Self, String> {
        if !hex(&inspection_revision) || fragment.index >= options.max_fragments {
            return Err("Choose a bounded fragment from a current digest inspection.".into());
        }
        let preview_context = hash(&(&selection, &options, &inspection_revision, &fragment))?;
        Ok(Self {
            source,
            selection,
            options,
            fragment,
            inspection_revision,
            preview_context,
        })
    }
    /// Accept a fresh native generic inspection, retaining its actual fragment
    /// contour; no Run-only adapter or imported duplex summary is involved.
    pub fn from_preview(
        source: SourcePin,
        molecule: &Molecule,
        selection: CatalogSelection,
        options: InspectionOptions,
        preview: &native::DigestPreview,
        fragment_index: usize,
    ) -> Result<Self, String> {
        source.validate_current(source.event(), source.epoch(), molecule, source.origin())?;
        let role = collagene_core::restriction_history::source_role(source.event())
            .map_err(|e| e.to_string())?;
        let expected = native::digest_inspection_revision(
            &role,
            molecule,
            &source.native_origin(),
            &selection,
            &options,
        )
        .map_err(|e| e.to_string())?;
        preview
            .validate_binding(&expected)
            .map_err(|e| e.to_string())?;
        if hash(&preview.selection)? != hash(&selection)? || !preview.inspection.products_determined
        {
            return Err(
                "The current native digest does not determine this catalog fragment.".into(),
            );
        }
        let fragment = preview
            .inspection
            .fragments
            .iter()
            .find(|fragment| fragment.index == fragment_index)
            .cloned()
            .ok_or("The selected native digest fragment is unavailable.")?;
        Self::capture(source, selection, options, expected, fragment)
    }
    pub fn source(&self) -> &SourcePin {
        &self.source
    }
    pub fn inspection_revision(&self) -> &str {
        &self.inspection_revision
    }
    pub fn validate_current(
        &self,
        event: &EventBinding,
        epoch: u64,
        source: &Molecule,
        origin: Option<&MaterialOrigin>,
        selection: &CatalogSelection,
        options: &InspectionOptions,
        inspection_revision: &str,
        fragment: &FragmentSummary,
    ) -> Result<(), String> {
        self.source.validate_current(event, epoch, source, origin)?;
        if hash(&(selection, options, inspection_revision, fragment))? != self.preview_context {
            return Err("Digest catalog, limits, inspection or fragment geometry changed; run the digest again.".into());
        }
        Ok(())
    }
    pub fn request(
        &self,
        workspace: &Workspace,
        name: &str,
    ) -> Result<PendingInstallation, String> {
        self.source.attest(workspace)?;
        PendingInstallation::new(
            workspace,
            &self.source,
            native::Request {
                source: self.source.event.clone(),
                origin: self.source.native_origin(),
                name: name.into(),
                save_primers: false,
                derivation: Derivation::CatalogDigest {
                    selection: self.selection.clone(),
                    options: self.options.clone(),
                    fragment_index: self.fragment.index,
                },
            },
        )
    }
}
/// A literal mismatch-PCR choice bound to its exact native source and policy.
#[derive(Debug, Clone)]
pub struct MismatchChoice {
    source: SourcePin,
    forward: Primer,
    reverse: Primer,
    options: collagene_core::mutagenesis::PcrOptions,
    selection: native::MismatchPcrSelection,
    preview_context: String,
}
impl MismatchChoice {
    pub fn from_preview(
        source: SourcePin,
        molecule: &Molecule,
        forward: Primer,
        reverse: Primer,
        options: collagene_core::mutagenesis::PcrOptions,
        preview: &native::MismatchPreview,
        selected_index: usize,
    ) -> Result<Self, String> {
        source.validate_current(source.event(), source.epoch(), molecule, source.origin())?;
        let context = native::mismatch_pcr_preview_revision(
            molecule,
            &source.native_origin(),
            &forward,
            &reverse,
            &options,
        )
        .map_err(|e| e.to_string())?;
        preview
            .validate_binding(&context)
            .map_err(|e| e.to_string())?;
        if preview.source_revision != molecule.revision()
            || preview.preview.report.source_revision != molecule.revision()
        {
            return Err("Mismatch-PCR preview names another source revision.".into());
        }
        let geometry = preview
            .preview
            .report
            .products
            .get(selected_index)
            .ok_or("Choose a reported mismatch-PCR geometry.")?
            .clone();
        let document = preview
            .preview
            .documents
            .get(selected_index)
            .ok_or("Mismatch-PCR geometry has no corresponding reference document.")?;
        if selected_index >= 10_000 || document.sequence != geometry.sequence {
            return Err("Mismatch-PCR preview document and geometry disagree.".into());
        }
        Ok(Self {
            source,
            forward,
            reverse,
            options,
            selection: native::MismatchPcrSelection {
                selected_index,
                expected_preview_revision: document.revision(),
                geometry,
            },
            preview_context: context,
        })
    }
    pub fn validate_current(
        &self,
        event: &EventBinding,
        epoch: u64,
        molecule: &Molecule,
        origin: Option<&MaterialOrigin>,
        forward: &Primer,
        reverse: &Primer,
        options: &collagene_core::mutagenesis::PcrOptions,
    ) -> Result<(), String> {
        self.source
            .validate_current(event, epoch, molecule, origin)?;
        if native::mismatch_pcr_preview_revision(
            molecule,
            &self.source.native_origin(),
            forward,
            reverse,
            options,
        )
        .map_err(|e| e.to_string())?
            != self.preview_context
        {
            return Err(
                "Mismatch-PCR primers or policy changed; review the current preview again.".into(),
            );
        }
        Ok(())
    }
    pub fn request(
        &self,
        workspace: &Workspace,
        name: &str,
        save_primers: bool,
    ) -> Result<PendingInstallation, String> {
        self.source.attest(workspace)?;
        PendingInstallation::new(
            workspace,
            &self.source,
            native::Request {
                source: self.source.event.clone(),
                origin: self.source.native_origin(),
                derivation: Derivation::MismatchPcr {
                    forward: self.forward.clone(),
                    reverse: self.reverse.clone(),
                    options: self.options.clone(),
                    selection: self.selection.clone(),
                },
                name: name.into(),
                save_primers,
            },
        )
    }
}
/// Repair acts on freshly reconstructed saved material, never imported end summaries.
#[derive(Debug, Clone)]
pub struct RepairChoice {
    source: SourcePin,
    method: collagene_core::chemistry::EndRepairMethod,
    preview_context: String,
}
impl RepairChoice {
    pub fn from_preview(
        source: SourcePin,
        molecule: &Molecule,
        method: collagene_core::chemistry::EndRepairMethod,
        preview: &native::RepairPreview,
    ) -> Result<Self, String> {
        source.validate_current(source.event(), source.epoch(), molecule, source.origin())?;
        let context = native::repair_preview_revision(molecule, &source.native_origin(), method)
            .map_err(|e| e.to_string())?;
        preview
            .validate_binding(&context)
            .map_err(|e| e.to_string())?;
        if preview.source_revision != molecule.revision() {
            return Err("End-repair preview names another source revision.".into());
        }
        Ok(Self {
            source,
            method,
            preview_context: context,
        })
    }
    pub fn validate_current(
        &self,
        event: &EventBinding,
        epoch: u64,
        molecule: &Molecule,
        origin: Option<&MaterialOrigin>,
        method: collagene_core::chemistry::EndRepairMethod,
    ) -> Result<(), String> {
        self.source
            .validate_current(event, epoch, molecule, origin)?;
        if native::repair_preview_revision(molecule, &self.source.native_origin(), method)
            .map_err(|e| e.to_string())?
            != self.preview_context
        {
            return Err(
                "The selected end-repair method changed; review the current chemistry again."
                    .into(),
            );
        }
        Ok(())
    }
    pub fn request(
        &self,
        workspace: &Workspace,
        name: &str,
    ) -> Result<PendingInstallation, String> {
        self.source.attest(workspace)?;
        PendingInstallation::new(
            workspace,
            &self.source,
            native::Request {
                source: self.source.event.clone(),
                origin: self.source.native_origin(),
                derivation: Derivation::EndRepair {
                    method: self.method,
                },
                name: name.into(),
                save_primers: false,
            },
        )
    }
}

/// Ready for the worker, not an installed product. The UI must retain its choice
/// and validate current context before dispatch, persistence and visible apply.
#[derive(Debug, Clone)]
pub struct PendingInstallation {
    request: native::Request,
    context_revision: String,
    workspace_fingerprint: String,
    source: SourcePin,
}
impl PendingInstallation {
    fn new(
        workspace: &Workspace,
        source: &SourcePin,
        request: native::Request,
    ) -> Result<Self, String> {
        native::preflight(workspace, &request).map_err(|e| e.to_string())?;
        Ok(Self {
            context_revision: native::revision(workspace, &request).map_err(|e| e.to_string())?,
            workspace_fingerprint: collagene_core::run_lineage::workspace_fingerprint(workspace)
                .map_err(|e| e.to_string())?,
            request,
            source: source.clone(),
        })
    }
    pub fn request(&self) -> &native::Request {
        &self.request
    }
    pub fn context_revision(&self) -> &str {
        &self.context_revision
    }
    pub fn validate_workspace(&self, workspace: &Workspace, epoch: u64) -> Result<(), String> {
        if epoch != self.source.epoch
            || collagene_core::run_lineage::workspace_fingerprint(workspace)
                .map_err(|e| e.to_string())?
                != self.workspace_fingerprint
        {
            return Err("The workspace changed while reusable material was being prepared.".into());
        }
        Ok(())
    }
    /// Validates the worker's archive receipt only. Caller still must atomically
    /// commit current+previous recovery frames and then visibly apply the archive.
    pub fn validate_received(
        &self,
        installed: &native::Installed,
        decoded: &Workspace,
    ) -> Result<(), String> {
        installed
            .validate_binding(&self.context_revision)
            .map_err(|e| e.to_string())?;
        if installed.receipt.source != self.source.event {
            return Err("Derivation receipt names another source event.".into());
        }
        installed
            .receipt
            .validate_archive(decoded)
            .map_err(|e| e.to_string())?;
        let actual = collagene_core::workspace::load_workspace(&installed.archive)
            .map_err(|e| e.to_string())?;
        if hash(&actual)? != hash(decoded)? {
            return Err(
                "Decoded workspace differs from the exact received derivation archive.".into(),
            );
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::{Topology, catalog_digest, restriction_history, verification};
    fn event(workspace: &Workspace) -> EventBinding {
        let doc = &workspace.documents[0];
        EventBinding {
            document_id: doc.id.clone(),
            event_id: doc.head_event.clone().unwrap(),
            revision: doc.head.clone(),
        }
    }
    fn pcr() -> (
        Workspace,
        Molecule,
        EventBinding,
        Primer,
        Primer,
        Arc<PcrDocument>,
        PcrChoice,
    ) {
        let source =
            Molecule::dna("Original template", "AGTCCGATGCTAGACT", Topology::Linear).unwrap();
        let workspace =
            Workspace::from_documents("Private workspace", vec![source.clone()], 0).unwrap();
        let source_event = event(&workspace);
        let forward = Primer::new("Forward", "AGTCCG", "").unwrap();
        let reverse = Primer::new("Reverse", "AGTCTA", "").unwrap();
        let previews = verification::simulate_pcr_documents(&source, &forward, &reverse).unwrap();
        assert_eq!(previews.len(), 1);
        let preview = Arc::new(previews.into_iter().next().unwrap());
        assert_eq!(preview.molecule.sequence.as_str(), "AGTCCGATGCTAGACT");
        let pin = SourcePin::capture(source_event.clone(), 7, &source, None).unwrap();
        let choice =
            PcrChoice::capture(pin, forward.clone(), reverse.clone(), 0, preview.clone()).unwrap();
        (
            workspace,
            source,
            source_event,
            forward,
            reverse,
            preview,
            choice,
        )
    }
    #[test]
    fn pcr_choice_binds_event_epoch_full_oligos_geometry_and_preview_not_equal_dna() {
        let (workspace, source, event, forward, reverse, preview, choice) = pcr();
        choice
            .validate_current(&event, 7, &source, None, &forward, &reverse, 0, &preview)
            .unwrap();
        let mut alias = event.clone();
        alias.document_id = "different-equal-DNA-document".into();
        assert!(
            choice
                .validate_current(&alias, 7, &source, None, &forward, &reverse, 0, &preview)
                .is_err()
        );
        assert!(
            choice
                .validate_current(&event, 8, &source, None, &forward, &reverse, 0, &preview)
                .is_err()
        );
        let tailed = Primer::new("Forward", "AGTCCG", "AAA").unwrap();
        assert!(
            choice
                .validate_current(&event, 7, &source, None, &tailed, &reverse, 0, &preview)
                .is_err()
        );
        assert!(
            choice
                .validate_current(&event, 7, &source, None, &forward, &reverse, 1, &preview)
                .is_err()
        );
        let mut changed = preview.as_ref().clone();
        changed.geometry.forward_start += 1;
        assert!(
            choice
                .validate_current(&event, 7, &source, None, &forward, &reverse, 0, &changed)
                .is_err()
        );
        let mut changed = preview.as_ref().clone();
        changed.molecule.name = "another preview".into();
        assert!(
            choice
                .validate_current(&event, 7, &source, None, &forward, &reverse, 0, &changed)
                .is_err()
        );
        let pending = choice.request(&workspace, "Reusable PCR").unwrap();
        assert!(matches!(
            pending.request().derivation,
            Derivation::ExactPcr { .. }
        ));
        pending.validate_workspace(&workspace, 7).unwrap();
        assert!(pending.validate_workspace(&workspace, 8).is_err());
        let mut changed = workspace;
        changed.name = "Different workspace context".into();
        assert!(pending.validate_workspace(&changed, 7).is_err());
    }
    #[test]
    fn received_native_pcr_archive_is_validated_but_no_commit_or_reference_promotion_occurs() {
        let (workspace, _, _, _, _, preview, choice) = pcr();
        let before = hash(&workspace).unwrap();
        let pending = choice
            .request_with_primers(&workspace, "Reusable PCR", true)
            .unwrap();
        let installed = native::install(&workspace, pending.request()).unwrap();
        assert!(installed.receipt.primer_attachment.is_some());
        let decoded = collagene_core::workspace::load_workspace(&installed.archive).unwrap();
        pending.validate_received(&installed, &decoded).unwrap();
        assert_eq!(hash(&workspace).unwrap(), before);
        assert_eq!(choice.preview_revision(), preview.molecule.revision());
        assert_ne!(
            installed.receipt.product.event_id,
            choice.source().event().event_id
        );
        let mut altered = decoded;
        altered.name = "Not the received archive".into();
        assert!(pending.validate_received(&installed, &altered).is_err());
        let mut wrong = installed;
        wrong.receipt.context_revision = "0".repeat(64);
        assert!(pending.validate_received(&wrong, &altered).is_err());
    }
    #[test]
    fn installed_derivation_origin_is_explicit_and_cannot_be_cleared_or_rebound_to_equal_dna() {
        let (workspace, _, _, _, _, _, choice) = pcr();
        let pending = choice.request(&workspace, "Reusable PCR").unwrap();
        let installed = native::install(&workspace, pending.request()).unwrap();
        let decoded = collagene_core::workspace::load_workspace(&installed.archive).unwrap();
        let product = installed.receipt.product.clone();
        let source = &decoded.revisions[&product.revision];
        let origin = Arc::new(MaterialOrigin::Derived(
            native::capsule(&decoded, &product).unwrap(),
        ));
        let pin = SourcePin::capture(product.clone(), 10, source, Some(origin.clone())).unwrap();
        pin.validate_current(&product, 10, source, Some(&origin))
            .unwrap();
        assert!(pin.validate_current(&product, 10, source, None).is_err());
        let mut alias = product.clone();
        alias.event_id = "equal-DNA-other-event".into();
        assert!(SourcePin::capture(alias, 10, source, Some(origin)).is_err());
        assert_eq!(pin.event(), &product);
    }
    #[test]
    fn catalog_fragment_retains_exact_geometry_limits_and_catalog_without_a_target() {
        let source = Molecule::dna("Original circle", "AAAGAATTCTTT", Topology::Circular).unwrap();
        let workspace =
            Workspace::from_documents("Digest context", vec![source.clone()], 0).unwrap();
        let event = event(&workspace);
        let pin = SourcePin::capture(event.clone(), 9, &source, None).unwrap();
        let catalog = catalog_digest::builtin_catalog().unwrap();
        let enzyme = catalog
            .edition()
            .enzymes
            .iter()
            .find(|enzyme| enzyme.name == "EcoRI")
            .unwrap();
        let selection = catalog_digest::select(&catalog, &[enzyme.id.clone()]).unwrap();
        let options = InspectionOptions::default();
        let role = restriction_history::source_role(&event).unwrap();
        let inspection =
            catalog_digest::inspect(&catalog, &selection, &role, &source, &options).unwrap();
        let fragment = inspection.inspection.fragments[0].clone();
        let choice = DigestChoice::capture(
            pin,
            selection.clone(),
            options.clone(),
            inspection.batch_revision.clone(),
            fragment.clone(),
        )
        .unwrap();
        choice
            .validate_current(
                &event,
                9,
                &source,
                None,
                &selection,
                &options,
                &inspection.batch_revision,
                &fragment,
            )
            .unwrap();
        let mut changed = options.clone();
        changed.max_sites -= 1;
        assert!(
            choice
                .validate_current(
                    &event,
                    9,
                    &source,
                    None,
                    &selection,
                    &changed,
                    &inspection.batch_revision,
                    &fragment
                )
                .is_err()
        );
        let mut changed = fragment;
        changed.top_start += 1;
        assert!(
            choice
                .validate_current(
                    &event,
                    9,
                    &source,
                    None,
                    &selection,
                    &options,
                    &inspection.batch_revision,
                    &changed
                )
                .is_err()
        );
        let pending = choice
            .request(&workspace, "Reusable digest fragment")
            .unwrap();
        assert!(matches!(
            pending.request().derivation,
            Derivation::CatalogDigest { .. }
        ));
        let installed = native::install(&workspace, pending.request()).unwrap();
        let decoded = collagene_core::workspace::load_workspace(&installed.archive).unwrap();
        pending.validate_received(&installed, &decoded).unwrap();
    }
    #[test]
    fn repeated_digest_uses_derived_origin_and_actual_generic_inspection() {
        let source = Molecule::dna(
            "Two-site original circle",
            "AAAGAATTCCCCGGGTTT",
            Topology::Circular,
        )
        .unwrap();
        let workspace =
            Workspace::from_documents("Repeated digest", vec![source.clone()], 0).unwrap();
        let original_event = event(&workspace);
        let catalog = catalog_digest::builtin_catalog().unwrap();
        let select = |name: &str| {
            let enzyme = catalog
                .edition()
                .enzymes
                .iter()
                .find(|enzyme| enzyme.name == name)
                .unwrap();
            catalog_digest::select(&catalog, &[enzyme.id.clone()]).unwrap()
        };
        let options = InspectionOptions::default();
        let role = restriction_history::source_role(&original_event).unwrap();
        let first_selection = select("EcoRI");
        let preview = native::inspect_digest(
            &role,
            &source,
            &SourceOrigin::OriginalReference,
            &first_selection,
            &options,
        )
        .unwrap();
        let first = DigestChoice::from_preview(
            SourcePin::capture(original_event, 1, &source, None).unwrap(),
            &source,
            first_selection,
            options.clone(),
            &preview,
            0,
        )
        .unwrap();
        let request = first.request(&workspace, "Opened EcoRI circle").unwrap();
        let installed = native::install(&workspace, request.request()).unwrap();
        let next = collagene_core::workspace::load_workspace(&installed.archive).unwrap();
        request.validate_received(&installed, &next).unwrap();
        let product = installed.receipt.product.clone();
        let material = &next.revisions[&product.revision];
        let origin = Arc::new(MaterialOrigin::Derived(
            native::capsule(&next, &product).unwrap(),
        ));
        let role = restriction_history::source_role(&product).unwrap();
        let second_selection = select("SmaI");
        let preview = native::inspect_digest(
            &role,
            material,
            &SourceOrigin::Saved(origin.as_ref().clone()),
            &second_selection,
            &options,
        )
        .unwrap();
        assert_eq!(preview.inspection.fragments.len(), 2);
        let second = DigestChoice::from_preview(
            SourcePin::capture(product, 2, material, Some(origin)).unwrap(),
            material,
            second_selection,
            options,
            &preview,
            0,
        )
        .unwrap();
        let request = second.request(&next, "Sequential SmaI fragment").unwrap();
        assert!(matches!(request.request().origin, SourceOrigin::Saved(_)));
        let installed = native::install(&next, request.request()).unwrap();
        let final_workspace =
            collagene_core::workspace::load_workspace(&installed.archive).unwrap();
        request
            .validate_received(&installed, &final_workspace)
            .unwrap();
        assert_ne!(installed.receipt.source, installed.receipt.product);
    }
}

#[cfg(test)]
mod mismatch_repair_tests {
    use super::*;
    use collagene_core::{Topology, chemistry::EndRepairMethod, mutagenesis::PcrOptions};
    fn binding(frame: &Workspace) -> EventBinding {
        let doc = frame
            .documents
            .iter()
            .find(|d| d.id == frame.active_document)
            .unwrap();
        EventBinding {
            document_id: doc.id.clone(),
            event_id: doc.head_event.clone().unwrap(),
            revision: doc.head.clone(),
        }
    }
    #[test]
    fn literal_mismatch_product_retains_oligos_and_policy_without_an_independent_target() {
        let source =
            Molecule::dna("Original template", "AGTCCGATGCTAGACT", Topology::Linear).unwrap();
        let frame = Workspace::from_documents("Literal mismatch", vec![source.clone()], 0).unwrap();
        let event = binding(&frame);
        let forward = Primer::new("Changed forward", "AATCCG", "").unwrap();
        let reverse = Primer::new("Reverse", "AGTCTA", "").unwrap();
        let options = PcrOptions {
            binding: collagene_core::primers::BindingOptions {
                max_mismatches: 1,
                min_three_prime_exact: 2,
                ..Default::default()
            },
            ..Default::default()
        };
        let preview = native::preview_mismatch_pcr(
            &source,
            &SourceOrigin::OriginalReference,
            &forward,
            &reverse,
            &options,
        )
        .unwrap();
        assert!(preview.reusable());
        assert_eq!(
            preview.preview.documents[0].sequence.as_str(),
            "AATCCGATGCTAGACT"
        );
        let pin = SourcePin::capture(event.clone(), 11, &source, None).unwrap();
        let choice = MismatchChoice::from_preview(
            pin.clone(),
            &source,
            forward.clone(),
            reverse.clone(),
            options.clone(),
            &preview,
            0,
        )
        .unwrap();
        let pending = choice.request(&frame, "Explicit mutant", true).unwrap();
        assert!(
            matches!(&pending.request().derivation,Derivation::MismatchPcr{selection,..} if selection.selected_index==0 && selection.geometry.sequence.as_str()=="AATCCGATGCTAGACT")
        );
        let before = serde_json::to_vec(&frame).unwrap();
        let installed = native::prepare_install(&frame, pending.request()).unwrap();
        assert!(installed.receipt.primer_attachment.is_some());
        assert_eq!(
            installed.workspace.revisions[&installed.receipt.product.revision]
                .sequence
                .as_str(),
            "AATCCGATGCTAGACT"
        );
        assert_eq!(serde_json::to_vec(&frame).unwrap(), before);
        assert!(
            choice
                .validate_current(&event, 12, &source, None, &forward, &reverse, &options)
                .is_err()
        );
        let mut changed = event.clone();
        changed.event_id = "equal-DNA-new-event".into();
        assert!(
            choice
                .validate_current(&changed, 11, &source, None, &forward, &reverse, &options)
                .is_err()
        );
        let mut changed_options = options.clone();
        changed_options.binding.min_three_prime_exact = 1;
        assert!(
            choice
                .validate_current(
                    &event,
                    11,
                    &source,
                    None,
                    &forward,
                    &reverse,
                    &changed_options
                )
                .is_err()
        );
        assert!(
            MismatchChoice::from_preview(pin, &source, forward, reverse, options, &preview, 1)
                .is_err()
        );
    }
    #[test]
    fn saved_sticky_fragment_is_repaired_directly_with_exact_origin_and_method_pins() {
        let source = Molecule::dna("EcoRI parent", "AAAGAATTCTTT", Topology::Linear).unwrap();
        let frame =
            Workspace::from_documents("Repair saved fragment", vec![source.clone()], 0).unwrap();
        let catalog = collagene_core::catalog_digest::builtin_catalog().unwrap();
        let enzyme = catalog
            .edition()
            .enzymes
            .iter()
            .find(|e| e.name == "EcoRI")
            .unwrap();
        let selection =
            collagene_core::catalog_digest::select(&catalog, &[enzyme.id.clone()]).unwrap();
        let cut = native::prepare_install(
            &frame,
            &native::Request {
                source: binding(&frame),
                origin: SourceOrigin::OriginalReference,
                derivation: Derivation::CatalogDigest {
                    selection,
                    options: Default::default(),
                    fragment_index: 0,
                },
                name: "Sticky fragment".into(),
                save_primers: false,
            },
        )
        .unwrap();
        let event = cut.receipt.product.clone();
        let molecule = cut.workspace.revisions[&event.revision].clone();
        assert_eq!(molecule.sequence.as_str(), "AAAG");
        let origin = Arc::new(MaterialOrigin::Derived(
            native::capsule(&cut.workspace, &event).unwrap(),
        ));
        let method = EndRepairMethod::FillFivePrimeOverhangs;
        let preview = native::preview_repair(
            &molecule,
            &SourceOrigin::Saved(origin.as_ref().clone()),
            method,
        )
        .unwrap();
        assert_eq!(preview.molecule.sequence.as_str(), "AAAGAATT");
        assert_eq!(preview.source_duplex.top_length, 4);
        assert_eq!(preview.source_duplex.bottom_length, 8);
        assert_eq!(preview.repair.duplex.top.len(), 8);
        assert_eq!(preview.repair.duplex.bottom.len(), 8);
        let pin = SourcePin::capture(event.clone(), 4, &molecule, Some(origin.clone())).unwrap();
        let choice = RepairChoice::from_preview(pin.clone(), &molecule, method, &preview).unwrap();
        let pending = choice
            .request(&cut.workspace, "Filled saved fragment")
            .unwrap();
        assert!(matches!(
            pending.request().derivation,
            Derivation::EndRepair {
                method: EndRepairMethod::FillFivePrimeOverhangs
            }
        ));
        let before = serde_json::to_vec(&cut.workspace).unwrap();
        let installed = native::prepare_install(&cut.workspace, pending.request()).unwrap();
        assert_eq!(
            installed.workspace.revisions[&installed.receipt.product.revision]
                .sequence
                .as_str(),
            "AAAGAATT"
        );
        assert_eq!(serde_json::to_vec(&cut.workspace).unwrap(), before);
        assert!(
            choice
                .validate_current(&event, 4, &molecule, None, method)
                .is_err()
        );
        assert!(
            choice
                .validate_current(
                    &event,
                    4,
                    &molecule,
                    Some(origin.as_ref()),
                    EndRepairMethod::RemoveThreePrimeOverhangs
                )
                .is_err()
        );
        assert!(
            RepairChoice::from_preview(
                pin,
                &molecule,
                EndRepairMethod::RemoveThreePrimeOverhangs,
                &preview
            )
            .is_err()
        );
        let mut other = event;
        other.event_id = "another-retained-event".into();
        assert!(SourcePin::capture(other, 4, &molecule, Some(origin)).is_err());
    }
}
