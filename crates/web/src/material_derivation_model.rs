//! Pure source/preview pins for explicit reusable PCR/digest installation.
//! This module never edits a workspace, writes recovery storage or reports commit.
use collagene_core::{
    Molecule,
    analysis::Primer,
    material_derivation::{self as native, Derivation, PcrSelection, SourceOrigin},
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
        Ok(Self {
            source,
            forward,
            reverse,
            selected_index,
            preview,
            preview_revision,
        })
    }

    pub fn request(
        &self,
        workspace: &Workspace,
        name: &str,
        save_primers: bool,
    ) -> Result<native::Request, String> {
        self.source.attest(workspace)?;
        let request = native::Request {
            source: self.source.event.clone(),
            origin: self.source.native_origin(),
            name: name.into(),
            save_primers,
            derivation: Derivation::ExactPcr {
                template_strand: self.preview.geometry.template_strand,
                forward: self.forward.clone(),
                reverse: self.reverse.clone(),
                selection: PcrSelection {
                    selected_index: self.selected_index,
                    expected_preview_revision: self.preview_revision.clone(),
                    geometry: self.preview.geometry.clone(),
                },
            },
        };
        native::preflight(workspace, &request).map_err(|e| e.to_string())?;
        Ok(request)
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
        })
    }

    pub fn request(
        &self,
        workspace: &Workspace,
        name: &str,
        save_primers: bool,
    ) -> Result<native::Request, String> {
        self.source.attest(workspace)?;
        let request = native::Request {
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
        };
        native::preflight(workspace, &request).map_err(|e| e.to_string())?;
        Ok(request)
    }
}
/// Repair acts on freshly reconstructed saved material, never imported end summaries.
#[derive(Debug, Clone)]
pub struct RepairChoice {
    source: SourcePin,
    method: collagene_core::chemistry::EndRepairMethod,
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
        Ok(Self { source, method })
    }

    pub fn request(&self, workspace: &Workspace, name: &str) -> Result<native::Request, String> {
        self.source.attest(workspace)?;
        let request = native::Request {
            source: self.source.event.clone(),
            origin: self.source.native_origin(),
            derivation: Derivation::EndRepair {
                method: self.method,
            },
            name: name.into(),
            save_primers: false,
        };
        native::preflight(workspace, &request).map_err(|e| e.to_string())?;
        Ok(request)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::{
        Topology, catalog_digest, restriction_history, restriction_inspection::InspectionOptions,
        verification,
    };
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
        let previews = verification::simulate_pcr_documents(
            &source,
            &forward,
            &reverse,
            collagene_core::chemistry::Strand::Top,
        )
        .unwrap();
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
    fn pcr_source_pins_and_native_geometry_refuse_stale_or_equal_dna_aliases() {
        let (workspace, source, event, _, _, _, choice) = pcr();
        choice
            .source
            .validate_current(&event, 7, &source, None)
            .unwrap();
        let mut alias = event.clone();
        alias.document_id = "different-equal-DNA-document".into();
        assert!(
            choice
                .source
                .validate_current(&alias, 7, &source, None)
                .is_err()
        );
        assert!(
            choice
                .source
                .validate_current(&event, 8, &source, None)
                .is_err()
        );
        let pending = choice.request(&workspace, "Reusable PCR", false).unwrap();
        let mut changed = pending.clone();
        if let Derivation::ExactPcr { forward, .. } = &mut changed.derivation {
            *forward = Primer::new("Forward", "AGTCCG", "AAA").unwrap();
        }
        assert!(native::prepare_install(&workspace, &changed).is_err());
        let mut changed = pending.clone();
        if let Derivation::ExactPcr { selection, .. } = &mut changed.derivation {
            selection.selected_index = 1;
        }
        assert!(native::prepare_install(&workspace, &changed).is_err());
        let mut changed = pending;
        if let Derivation::ExactPcr { selection, .. } = &mut changed.derivation {
            selection.geometry.forward_start += 1;
        }
        assert!(native::prepare_install(&workspace, &changed).is_err());
    }
    #[test]
    fn received_native_pcr_archive_is_validated_but_no_commit_or_reference_promotion_occurs() {
        let (workspace, _, _, _, _, preview, choice) = pcr();
        let before = hash(&workspace).unwrap();
        let pending = choice.request(&workspace, "Reusable PCR", true).unwrap();
        let installed = native::install(&workspace, &pending).unwrap();
        assert!(installed.receipt.primer_attachment.is_some());
        let decoded = collagene_core::workspace::load_workspace(&installed.archive).unwrap();
        installed
            .validate_binding(&native::revision(&workspace, &pending).unwrap())
            .unwrap();
        installed.receipt.validate_archive(&decoded).unwrap();
        assert_eq!(hash(&workspace).unwrap(), before);
        assert_eq!(choice.preview_revision, preview.molecule.revision());
        assert_ne!(
            installed.receipt.product.event_id,
            choice.source.event().event_id
        );
        let mut wrong = installed;
        wrong.receipt.context_revision = "0".repeat(64);
        assert!(
            wrong
                .validate_binding(&native::revision(&workspace, &pending).unwrap())
                .is_err()
        );
    }
    #[test]
    fn installed_derivation_origin_is_explicit_and_cannot_be_cleared_or_rebound_to_equal_dna() {
        let (workspace, _, _, _, _, _, choice) = pcr();
        let pending = choice.request(&workspace, "Reusable PCR", false).unwrap();
        let installed = native::install(&workspace, &pending).unwrap();
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
    fn catalog_preview_binds_current_limits_and_literal_cut_geometry_without_a_target() {
        let source = Molecule::dna("Original circle", "AAAGAATTCTTT", Topology::Circular).unwrap();
        let workspace = Workspace::from_documents("Digest", vec![source.clone()], 0).unwrap();
        let event = event(&workspace);
        let role = restriction_history::source_role(&event).unwrap();
        let catalog = catalog_digest::builtin_catalog().unwrap();
        let selection = catalog_digest::select(&catalog, &["ecori".into()]).unwrap();
        let options = collagene_core::restriction_inspection::InspectionOptions::default();
        let preview = native::inspect_digest(
            &role,
            &source,
            &SourceOrigin::OriginalReference,
            &selection,
            &options,
        )
        .unwrap();
        let expected = native::digest_inspection_revision(
            &role,
            &source,
            &SourceOrigin::OriginalReference,
            &selection,
            &options,
        )
        .unwrap();
        preview.validate_binding(&expected).unwrap();
        let fragment = &preview.inspection.fragments[0];
        assert_eq!(fragment.top_start, 4);
        assert_eq!(fragment.bottom_start, 8);
        let mut changed = options.clone();
        changed.max_sites -= 1;
        let changed_context = native::digest_inspection_revision(
            &role,
            &source,
            &SourceOrigin::OriginalReference,
            &selection,
            &changed,
        )
        .unwrap();
        assert!(preview.validate_binding(&changed_context).is_err());
        let request = native::Request {
            source: event.clone(),
            origin: SourceOrigin::OriginalReference,
            name: "Reusable digest fragment".into(),
            save_primers: false,
            derivation: Derivation::CatalogDigest {
                selection,
                options,
                fragment_index: fragment.index,
            },
        };
        let installed = native::install(&workspace, &request).unwrap();
        installed
            .validate_binding(&native::revision(&workspace, &request).unwrap())
            .unwrap();
        let decoded = collagene_core::workspace::load_workspace(&installed.archive).unwrap();
        installed.receipt.validate_archive(&decoded).unwrap();
        assert_eq!(installed.receipt.source, event);
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
        preview
            .validate_binding(
                &native::digest_inspection_revision(
                    &role,
                    &source,
                    &SourceOrigin::OriginalReference,
                    &first_selection,
                    &options,
                )
                .unwrap(),
            )
            .unwrap();
        let request = native::Request {
            source: original_event,
            origin: SourceOrigin::OriginalReference,
            name: "Opened EcoRI circle".into(),
            save_primers: false,
            derivation: Derivation::CatalogDigest {
                selection: first_selection,
                options: options.clone(),
                fragment_index: preview.inspection.fragments[0].index,
            },
        };
        let installed = native::install(&workspace, &request).unwrap();
        let next = collagene_core::workspace::load_workspace(&installed.archive).unwrap();
        installed
            .validate_binding(&native::revision(&workspace, &request).unwrap())
            .unwrap();
        installed.receipt.validate_archive(&next).unwrap();
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
        preview
            .validate_binding(
                &native::digest_inspection_revision(
                    &role,
                    material,
                    &SourceOrigin::Saved(origin.as_ref().clone()),
                    &second_selection,
                    &options,
                )
                .unwrap(),
            )
            .unwrap();
        let request = native::Request {
            source: product,
            origin: SourceOrigin::Saved(origin.as_ref().clone()),
            name: "Sequential SmaI fragment".into(),
            save_primers: false,
            derivation: Derivation::CatalogDigest {
                selection: second_selection,
                options,
                fragment_index: preview.inspection.fragments[0].index,
            },
        };
        assert!(matches!(&request.origin, SourceOrigin::Saved(_)));
        let installed = native::install(&next, &request).unwrap();
        let final_workspace =
            collagene_core::workspace::load_workspace(&installed.archive).unwrap();
        installed
            .validate_binding(&native::revision(&next, &request).unwrap())
            .unwrap();
        installed
            .receipt
            .validate_archive(&final_workspace)
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
            matches!(&pending.derivation,Derivation::MismatchPcr{selection,..} if selection.selected_index==0 && selection.geometry.sequence.as_str()=="AATCCGATGCTAGACT")
        );
        let before = serde_json::to_vec(&frame).unwrap();
        let installed = native::prepare_install(&frame, &pending).unwrap();
        assert!(installed.receipt.primer_attachment.is_some());
        assert_eq!(
            installed.workspace.revisions[&installed.receipt.product.revision]
                .sequence
                .as_str(),
            "AATCCGATGCTAGACT"
        );
        assert_eq!(serde_json::to_vec(&frame).unwrap(), before);
        assert!(pin.validate_current(&event, 12, &source, None).is_err());
        let mut changed = event.clone();
        changed.event_id = "equal-DNA-new-event".into();
        assert!(pin.validate_current(&changed, 11, &source, None).is_err());
        let mut changed_options = options.clone();
        changed_options.binding.min_three_prime_exact = 1;
        assert!(
            MismatchChoice::from_preview(
                pin.clone(),
                &source,
                forward.clone(),
                reverse.clone(),
                changed_options,
                &preview,
                0
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
            pending.derivation,
            Derivation::EndRepair {
                method: EndRepairMethod::FillFivePrimeOverhangs
            }
        ));
        let before = serde_json::to_vec(&cut.workspace).unwrap();
        let installed = native::prepare_install(&cut.workspace, &pending).unwrap();
        assert_eq!(
            installed.workspace.revisions[&installed.receipt.product.revision]
                .sequence
                .as_str(),
            "AAAGAATT"
        );
        assert_eq!(serde_json::to_vec(&cut.workspace).unwrap(), before);
        assert!(pin.validate_current(&event, 4, &molecule, None).is_err());
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
