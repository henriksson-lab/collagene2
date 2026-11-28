//! Pure PCR-product digest selection. Source events and UI epochs stay with the
//! caller; no original contribution interval is inferred from an amplicon.
use crate::cloning::preparation::cloning_methods_model::MethodsState;
use crate::cloning::preparation::cloning_methods_model::OverlapPreparation;
use collagene_core::{
    cloning::{
        pcr::digest::{Preview, Request},
        restriction::inspection::InspectionOptions,
    },
    sequence::analysis::Primer,
};

/// Match the native origin digest using bounded streaming JSON. This hashes a
/// documentary capsule only; it never restores physical state or replays it.
fn material_revision(
    origin: &collagene_core::material::run::MaterialOrigin,
) -> Result<String, String> {
    origin.validate_documentary().map_err(|e| e.to_string())?;
    origin_revision(origin)
}
fn origin_revision(origin: &impl serde::Serialize) -> Result<String, String> {
    use std::io::Write;
    struct Digest {
        hash: blake3::Hasher,
        bytes: usize,
    }
    impl Write for Digest {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.bytes = self
                .bytes
                .checked_add(bytes.len())
                .filter(|n| *n <= 32 << 20)
                .ok_or_else(|| std::io::Error::other("bounded origin digest"))?;
            self.hash.update(bytes);
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut writer = Digest {
        hash: blake3::Hasher::new(),
        bytes: 0,
    };
    serde_json::to_writer(&mut writer, origin).map_err(|e| e.to_string())?;
    Ok(writer.hash.finalize().to_hex().to_string())
}

pub(crate) fn build_request(state: &MethodsState, row: usize) -> Result<Request, String> {
    if state
        .rows
        .get(row)
        .is_none_or(|r| r.preparation != OverlapPreparation::PcrCatalogDigest)
    {
        return Err("Choose exact PCR followed by a catalog digest.".into());
    }
    build_base(state, row)
}
fn build_base(state: &MethodsState, row: usize) -> Result<Request, String> {
    let roles = state.role_ids(state.rows.len())?;
    let input = state
        .rows
        .get(row)
        .ok_or("The preparation row is unavailable.")?;
    if !matches!(
        input.preparation,
        OverlapPreparation::PcrCatalogDigest | OverlapPreparation::MismatchPcrCatalogDigest
    ) {
        return Err(
            "Choose PCR followed by a pinned catalog digest before inspecting fragments.".into(),
        );
    }
    for name in [&input.forward_name, &input.reverse_name] {
        if name.trim().is_empty() || name.len() > 256 || name.chars().any(char::is_control) {
            return Err("Each PCR primer needs a nonblank name of at most256 bytes.".into());
        }
    }
    if [
        &input.forward,
        &input.reverse_primer,
        &input.forward_tail,
        &input.reverse_tail,
    ]
    .iter()
    .any(|text| text.len() > 65536)
    {
        return Err("PCR primer fields exceed64KiB.".into());
    }
    if input.fragment.len() > 20 {
        return Err("The digest-fragment field exceeds its bounded size.".into());
    }
    let forward = Primer::new(&input.forward_name, &input.forward, &input.forward_tail)
        .map_err(|error| error.to_string())?;
    let reverse = Primer::new(
        &input.reverse_name,
        &input.reverse_primer,
        &input.reverse_tail,
    )
    .map_err(|error| error.to_string())?;
    if forward.annealing.is_empty() || reverse.annealing.is_empty() {
        return Err(
            "Enter both nonempty annealing domains before inspecting the PCR product.".into(),
        );
    }
    for primer in [&forward, &reverse] {
        if primer.annealing.len() > 1000
            || primer.tail.len() > 1000
            || !primer
                .annealing
                .as_str()
                .bytes()
                .chain(primer.tail.bytes())
                .all(|base| b"ACGT".contains(&base))
        {
            return Err("PCR-product inspection requires concrete annealing domains/tails of at most1000 bases each.".into());
        }
    }
    let selection = input
        .catalog_selection
        .clone()
        .ok_or("Choose an exact pinned catalog enzyme selection first.")?;
    selection.validate().map_err(|error| error.to_string())?;
    Ok(Request {
        template_strand: if input.preparation == OverlapPreparation::MismatchPcrCatalogDigest {
            collagene_core::chemistry::Strand::Top
        } else {
            input.template_strand
        },
        source_role: roles[row].role.clone(),
        forward,
        reverse,
        selection,
        options: InspectionOptions::default(),
    })
}

fn same_primer(a: &Primer, b: &Primer) -> bool {
    a.name == b.name && a.annealing == b.annealing && a.tail == b.tail
}
fn same_request(a: &Request, b: &Request) -> bool {
    a.template_strand == b.template_strand
        && a.source_role == b.source_role
        && same_primer(&a.forward, &b.forward)
        && same_primer(&a.reverse, &b.reverse)
        && a.selection == b.selection
        && a.options == b.options
}

/// The host checks its exact original-source event and epoch before calling.
/// An explicit selection changes only the PCR-product fragment number.
pub(crate) fn apply(
    state: &MethodsState,
    row: usize,
    bound_request: &Request,
    preview: &Preview,
    expected_revision: &str,
    fragment: usize,
) -> Result<MethodsState, String> {
    let mut current = build_request(state, row)?;
    // Search limits belong to the explicitly captured inspection request, not
    // to editable primer rows. The native binding validator checks their bounds.
    current.options = bound_request.options.clone();
    if !same_request(&current, bound_request) || !same_request(&preview.request, bound_request) {
        return Err(
            "Primers, source role or catalog settings changed; inspect the PCR digest again."
                .into(),
        );
    }
    preview
        .validate_binding(expected_revision)
        .map_err(|error| error.to_string())?;
    if preview.pcr_geometries.len() != 1 || preview.amplicon.is_none() {
        return Err("A unique freshly reconstructed PCR product is required before choosing its digest fragment.".into());
    }
    let inspection = preview
        .inspection
        .as_ref()
        .ok_or("The unique PCR product has no physical digest inspection.")?;
    if !inspection.scan_complete || !inspection.products_determined {
        return Err("The PCR digest is incomplete or its physical fragments are undetermined; no selection was applied.".into());
    }
    if !inspection
        .fragments
        .iter()
        .any(|item| item.index == fragment)
    {
        return Err("The explicitly chosen PCR-product digest fragment is unavailable.".into());
    }
    let number = fragment
        .checked_add(1)
        .ok_or("Digest-fragment number overflow.")?;
    let mut next = state.clone();
    next.rows[row].fragment = number.to_string();
    Ok(next)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assembly::draft::assembly_role_ids::Role;
    use crate::cloning::preparation::cloning_methods_model::Mode;
    use crate::cloning::preparation::cloning_methods_model::Treatment;
    use collagene_core::cloning::{
        pcr::digest as pcr_digest_inspection, restriction::catalog as catalog_digest,
    };
    use collagene_model::{Molecule, Topology};

    fn select(name: &str) -> collagene_core::cloning::restriction::catalog::CatalogSelection {
        let catalog = catalog_digest::builtin_catalog().unwrap();
        let id = catalog
            .edition()
            .enzymes
            .iter()
            .find(|entry| entry.name == name)
            .unwrap()
            .id
            .clone();
        catalog_digest::select(&catalog, &[id]).unwrap()
    }
    fn fixture() -> (MethodsState, Molecule, Request, Preview) {
        let mut state = MethodsState::default();
        state.mode = Mode::Ligation;
        state.source_roles = Some(vec![
            Role {
                input_id: "backbone".into(),
                role: "retained-vector".into(),
            },
            Role {
                input_id: "original-template".into(),
                role: "insert-amplification".into(),
            },
        ]);
        let row = &mut state.rows[1];
        row.preparation = OverlapPreparation::PcrCatalogDigest;
        row.document = 17;
        row.interval_start = "4".into();
        row.interval_end = "17".into();
        row.reverse = true;
        row.treatment = Treatment::Phosphorylate;
        row.fragment = "29".into();
        row.right_fusion = "CGCG".into();
        row.forward_name = "Forward with EcoRI".into();
        row.forward = "ATGCCATG".into();
        row.forward_tail = "GGGAATTCCC".into();
        row.reverse_name = "Reverse with EcoRI".into();
        row.reverse_primer = "TACCGGTA".into();
        row.reverse_tail = "TTGAATTCGG".into();
        row.catalog_selection = Some(select("EcoRI"));
        let source = Molecule::dna(
            "Original template",
            "ATGCCATGTCAGTACCGGTA",
            Topology::Linear,
        )
        .unwrap();
        let request = build_request(&state, 1).unwrap();
        let preview = pcr_digest_inspection::inspect(&source, &request).unwrap();
        (state, source, request, preview)
    }

    #[test]
    fn derived_digest_application_keeps_origin_and_refuses_stale_or_removed_binding() {
        use collagene_core::{
            chemistry::Strand,
            documents::workspace::Workspace,
            material::{derivation as d, run::MaterialOrigin},
        };
        let (mut state, source, request, _) = fixture();
        let w = Workspace::from_documents("Native template derivation", vec![source.clone()], 0)
            .unwrap();
        let doc = &w.documents[0];
        let event = collagene_core::material::lineage::EventBinding {
            document_id: doc.id.clone(),
            event_id: doc.head_event.clone().unwrap(),
            revision: doc.head.clone(),
        };
        let forward = Primer::new("Native F", "ATGCCATG", "").unwrap();
        let reverse = Primer::new("Native R", "TACCGGTA", "").unwrap();
        let original = d::SourceOrigin::OriginalReference;
        let documents =
            d::pcr::preview_pcr(&source, &original, &forward, &reverse, Strand::Top).unwrap();
        assert_eq!(documents.len(), 1);
        let prepared = d::prepare_install(
            &w,
            &d::Request {
                sources: vec![d::BoundSource {
                    role: d::SourceRole::Template,
                    event: event,
                    origin: original.clone(),
                }],
                derivation: d::Derivation::ExactPcr {
                    selection: d::PcrSelection {
                        selected_index: 0,
                        expected_preview_revision: documents[0].molecule.revision(),
                        geometry: documents[0].geometry.clone(),
                    },
                    forward,
                    reverse,
                    template_strand: Strand::Top,
                },
                name: "Retained PCR template".into(),
                save_primers: false,
            },
        )
        .unwrap();
        let origin = MaterialOrigin::Derived(
            d::capsule(&prepared.workspace, &prepared.receipt.product).unwrap(),
        );
        let source = &prepared.workspace.revisions[&origin.event().revision];
        state.rows[1].run_origin = Some(
            crate::history::run_origin_model::PinnedOrigin::from_material(origin.clone()).unwrap(),
        );
        let report = pcr_digest_inspection::inspect_from_material(
            source,
            &origin,
            &request,
            &prepared.workspace,
        )
        .unwrap();
        let revision = report.context_revision.clone();
        let request = InspectionRequest::Exact(request);
        let report = InspectionPreview::Material(report);
        let selected = apply_inspection(&state, 1, &request, &report, &revision, 1).unwrap();
        assert_eq!(selected.rows[1].fragment, "2");
        assert!(selected.rows[1].run_origin == state.rows[1].run_origin);
        let mut changed = state.clone();
        changed.rows[1].forward_tail.push('A');
        assert!(apply_inspection(&changed, 1, &request, &report, &revision, 1).is_err());
        state.rows[1].run_origin = None;
        assert!(apply_inspection(&state, 1, &request, &report, &revision, 1).is_err());
    }

    #[test]
    fn request_retains_exact_role_names_domains_tails_and_catalog_pin() {
        let (state, _, request, _) = fixture();
        assert_eq!(request.source_role, "insert-amplification");
        assert_eq!(request.forward.name, "Forward with EcoRI");
        assert_eq!(request.forward.annealing.as_str(), "ATGCCATG");
        assert_eq!(request.forward.tail, "GGGAATTCCC");
        assert_eq!(request.reverse.name, "Reverse with EcoRI");
        assert_eq!(request.reverse.annealing.as_str(), "TACCGGTA");
        assert_eq!(request.reverse.tail, "TTGAATTCGG");
        assert_eq!(
            request.selection,
            state.rows[1].catalog_selection.clone().unwrap()
        );
        assert_eq!(request.options, InspectionOptions::default());
    }

    #[test]
    fn explicit_fragment_selection_changes_only_its_one_based_field() {
        let (state, source, request, preview) = fixture();
        let before = state.clone();
        assert_eq!(preview.pcr_geometries.len(), 1);
        let digest = preview.inspection.as_ref().unwrap();
        assert!(digest.scan_complete && digest.products_determined && digest.fragments.len() >= 3);
        let index = digest.fragments.last().unwrap().index;
        let expected = pcr_digest_inspection::inspection_revision(&source, &request).unwrap();
        let applied = apply(&state, 1, &request, &preview, &expected, index).unwrap();
        let mut wanted = state.clone();
        wanted.rows[1].fragment = (index + 1).to_string();
        assert!(applied == wanted);
        assert!(state == before);
        assert_eq!(applied.rows[1].interval_start, "4");
        assert_eq!(applied.rows[1].interval_end, "17");
        assert!(applied.rows[1].reverse && applied.rows[1].treatment == Treatment::Phosphorylate);
    }

    #[test]
    fn stale_primers_roles_catalog_and_wrong_revision_or_fragment_never_apply() {
        let (state, source, request, preview) = fixture();
        let revision = pcr_digest_inspection::inspection_revision(&source, &request).unwrap();
        let fragment = preview.inspection.as_ref().unwrap().fragments[0].index;
        for changed in [
            {
                let mut changed = state.clone();
                changed.rows[1].forward_tail.push('A');
                changed
            },
            {
                let mut changed = state.clone();
                changed.rows[1].reverse_name = "Renamed reverse".into();
                changed
            },
            {
                let mut changed = state.clone();
                changed.source_roles.as_mut().unwrap()[1].role = "different-role".into();
                changed
            },
            {
                let mut changed = state.clone();
                changed.rows[1].catalog_selection = Some(select("BamHI"));
                changed
            },
        ] {
            assert!(apply(&changed, 1, &request, &preview, &revision, fragment).is_err());
        }
        assert!(apply(&state, 1, &request, &preview, &"f".repeat(64), fragment).is_err());
        assert!(apply(&state, 1, &request, &preview, &revision, usize::MAX).is_err());
        let mut unrelated = preview.clone();
        unrelated.request.reverse.tail.push('A');
        assert!(apply(&state, 1, &request, &unrelated, &revision, fragment).is_err());
        let mut wrong = state.clone();
        wrong.rows[1].preparation = OverlapPreparation::CatalogDigest;
        assert!(build_request(&wrong, 1).is_err());
        assert!(build_request(&state, 12).is_err());
    }

    #[test]
    fn incomplete_digestion_and_ambiguous_pcr_do_not_choose_a_default() {
        let (state, source, request, _) = fixture();
        let mut limited = request.clone();
        limited.options.max_comparisons = 1;
        let partial = pcr_digest_inspection::inspect(&source, &limited).unwrap();
        assert!(!partial.inspection.as_ref().unwrap().scan_complete);
        let revision = pcr_digest_inspection::inspection_revision(&source, &limited).unwrap();
        assert!(apply(&state, 1, &limited, &partial, &revision, 0).is_err());
        let repeated = Molecule::dna(
            "Two matching template regions",
            "ATGCCATGTCAGTACCGGTAATGCCATGTCAGTACCGGTA",
            Topology::Linear,
        )
        .unwrap();
        let ambiguous = pcr_digest_inspection::inspect(&repeated, &request).unwrap();
        assert!(
            ambiguous.pcr_geometries.len() > 1
                && ambiguous.amplicon.is_none()
                && ambiguous.inspection.is_none()
        );
        let revision = pcr_digest_inspection::inspection_revision(&repeated, &request).unwrap();
        assert!(apply(&state, 1, &request, &ambiguous, &revision, 0).is_err());
    }
}

/// Distinct preparation reports retain their native request and source bindings.
#[derive(Clone, serde::Serialize)]
pub(crate) enum InspectionRequest {
    Exact(Request),
    Mismatch(collagene_core::cloning::pcr::mismatch_digest::Request),
    Extension(collagene_core::cloning::overlap::extension::digest::Request),
}
pub(crate) enum InspectionPreview {
    Exact(Preview),
    Material(collagene_core::cloning::pcr::digest::MaterialPreview),
    Mismatch(collagene_core::cloning::pcr::mismatch_digest::Preview),
    Extension(collagene_core::cloning::overlap::extension::digest::Preview),
}
pub(crate) fn build_inspection_request(
    state: &MethodsState,
    row: usize,
) -> Result<InspectionRequest, String> {
    if let Some(input) = state
        .rows
        .get(row)
        .filter(|input| input.preparation == OverlapPreparation::OverlapExtensionPcrCatalogDigest)
    {
        let roles = state.role_ids(state.rows.len())?;
        if input.fragment.len() > 20 {
            return Err("The digest-fragment field exceeds its bounded size.".into());
        }
        let preparation = input.overlap_extension_pcr.as_ref().ok_or("Enter all six overlap-extension oligos before inspecting the outer-PCR product digest.")?.request()?;
        let selection = input
            .catalog_selection
            .as_ref()
            .ok_or("Pin the outer-PCR product digest enzymes first.")?;
        selection.validate().map_err(|e| e.to_string())?;
        return Ok(InspectionRequest::Extension(
            collagene_core::cloning::overlap::extension::digest::Request {
                source_role: roles[row].role.clone(),
                preparation,
                selection: selection.clone(),
                options: InspectionOptions::default(),
            },
        ));
    }
    let base = build_base(state, row)?;
    if state.rows[row].preparation == OverlapPreparation::MismatchPcrCatalogDigest {
        Ok(InspectionRequest::Mismatch(
            collagene_core::cloning::pcr::mismatch_digest::Request {
                base,
                pcr_options: state.rows[row].mismatch_options()?,
            },
        ))
    } else {
        Ok(InspectionRequest::Exact(base))
    }
}
impl InspectionRequest {
    pub(crate) fn key(&self) -> Result<String, String> {
        serde_json::to_string(self).map_err(|e| e.to_string())
    }
    pub(crate) fn revision(
        &self,
        source: &collagene_model::Molecule,
        origin: Option<&collagene_core::material::run::MaterialOrigin>,
        workspace: Option<&collagene_core::documents::workspace::Workspace>,
    ) -> Result<String, collagene_model::Error> {
        use collagene_model::Error;
        match self {
            Self::Exact(r) => match origin {
                Some(o) => collagene_core::cloning::pcr::digest::inspection_revision_from_material(
                    source,
                    o,
                    r,
                    workspace.ok_or_else(|| {
                        Error::Unsupported(
                            "The exact owning workspace is required for saved-template inspection."
                                .into(),
                        )
                    })?,
                ),
                None => collagene_core::cloning::pcr::digest::inspection_revision(source, r),
            },
            Self::Mismatch(r) => {
                collagene_core::cloning::pcr::mismatch_digest::inspection_revision(
                    source, origin, r, workspace,
                )
            }
            Self::Extension(r) => {
                collagene_core::cloning::overlap::extension::digest::inspection_revision(
                    source, origin, r, workspace,
                )
            }
        }
    }
    pub(crate) fn into_job(
        self,
        molecule: collagene_model::Molecule,
        origin: Option<collagene_core::material::run::MaterialOrigin>,
        workspace: Option<collagene_core::documents::workspace::Workspace>,
    ) -> Result<collagene_core::protocol::schema::AnalysisJob, String> {
        use collagene_core::protocol::schema::AnalysisJob;
        Ok(match self {
            Self::Exact(request) => match origin {
                Some(origin) => AnalysisJob::InspectMaterialPcrDigest {
                    molecule,
                    origin,
                    request,
                    workspace: workspace.ok_or("The owning workspace is unavailable.")?,
                },
                None => AnalysisJob::InspectPcrDigest { molecule, request },
            },
            Self::Mismatch(request) => AnalysisJob::InspectMismatchPcrDigest {
                molecule,
                origin,
                request,
                workspace,
            },
            Self::Extension(request) => AnalysisJob::InspectOverlapExtensionDigest {
                source: molecule,
                origin,
                request,
                workspace,
            },
        })
    }
}

impl InspectionPreview {
    pub(crate) fn context_revision(&self) -> &str {
        match self {
            Self::Material(p) => &p.context_revision,
            Self::Exact(p) => &p.context_revision,
            Self::Mismatch(p) => &p.context_revision,
            Self::Extension(p) => &p.context_revision,
        }
    }
    pub(crate) fn source_revision(&self) -> &str {
        match self {
            Self::Material(p) => &p.preview.source_revision,
            Self::Exact(p) => &p.source_revision,
            Self::Mismatch(p) => &p.source_revision,
            Self::Extension(p) => &p.source_revision,
        }
    }
    pub(crate) fn request(&self) -> InspectionRequest {
        match self {
            Self::Material(p) => InspectionRequest::Exact(p.preview.request.clone()),
            Self::Exact(p) => InspectionRequest::Exact(p.request.clone()),
            Self::Mismatch(p) => InspectionRequest::Mismatch(p.request.clone()),
            Self::Extension(p) => InspectionRequest::Extension(p.request.clone()),
        }
    }
    pub(crate) fn validate_binding(&self, r: &str) -> Result<(), collagene_model::Error> {
        match self {
            Self::Material(p) => p.validate_binding(r),
            Self::Exact(p) => p.validate_binding(r),
            Self::Mismatch(p) => p.validate_binding(r),
            Self::Extension(p) => p.validate_binding(r),
        }
    }
    pub(crate) fn amplicon(&self) -> Option<&collagene_model::Molecule> {
        match self {
            Self::Material(p) => p.preview.amplicon.as_ref(),
            Self::Exact(p) => p.amplicon.as_ref(),
            Self::Mismatch(p) => p.amplicon.as_ref(),
            Self::Extension(p) => p.amplicon.as_ref(),
        }
    }
    pub(crate) fn inspection(
        &self,
    ) -> Option<&collagene_core::cloning::restriction::inspection::RestrictionInspection> {
        match self {
            Self::Material(p) => p.preview.inspection.as_ref(),
            Self::Exact(p) => p.inspection.as_ref(),
            Self::Mismatch(p) => p.inspection.as_ref(),
            Self::Extension(p) => p.inspection.as_ref(),
        }
    }
    pub(crate) fn findings(&self) -> &[String] {
        match self {
            Self::Material(p) => &p.preview.findings,
            Self::Exact(p) => &p.findings,
            Self::Mismatch(p) => &p.findings,
            Self::Extension(p) => &p.findings,
        }
    }
    pub(crate) fn geometry_count(&self) -> usize {
        match self {
            Self::Material(p) => p.preview.pcr_geometries.len(),
            Self::Exact(p) => p.pcr_geometries.len(),
            Self::Mismatch(p) => p.pcr.report.products.len(),
            Self::Extension(p) => usize::from(p.amplicon.is_some()),
        }
    }
    pub(crate) fn geometry_page(
        &self,
        start: usize,
    ) -> Vec<collagene_core::cloning::pcr::digest::GeometrySummary> {
        match self {
            Self::Material(p) => p
                .preview
                .pcr_geometries
                .iter()
                .skip(start)
                .take(20)
                .cloned()
                .collect(),
            Self::Extension(_) => Vec::new(),
            Self::Exact(p) => p
                .pcr_geometries
                .iter()
                .skip(start)
                .take(20)
                .cloned()
                .collect(),
            Self::Mismatch(p) => p
                .pcr
                .report
                .products
                .iter()
                .enumerate()
                .skip(start)
                .take(20)
                .map(
                    |(index, g)| collagene_core::cloning::pcr::digest::GeometrySummary {
                        template_strand: collagene_core::chemistry::Strand::Top,
                        index,
                        forward_start: g.forward_start,
                        reverse_start: g.reverse_start,
                        product_length: g.sequence.len(),
                    },
                )
                .collect(),
        }
    }
    pub(crate) fn incomplete(&self) -> bool {
        match self {
            Self::Material(p) => p.origin_incomplete,
            Self::Exact(_) => false,
            Self::Mismatch(p) => p.incomplete,
            Self::Extension(p) => p.incomplete,
        }
    }
    pub(crate) fn selectable(&self) -> bool {
        if let Self::Extension(p) = self {
            return p.selectable();
        }
        self.geometry_count() == 1
            && self.amplicon().is_some()
            && self
                .inspection()
                .is_some_and(|i| i.scan_complete && i.products_determined)
    }
}
/// Caller retains the exact original source event, saved origin and UI epoch.
pub(crate) fn apply_inspection(
    state: &MethodsState,
    row: usize,
    bound: &InspectionRequest,
    preview: &InspectionPreview,
    revision: &str,
    fragment: usize,
) -> Result<MethodsState, String> {
    match (bound, preview) {
        (InspectionRequest::Exact(r), InspectionPreview::Material(p)) => {
            p.validate_binding(revision).map_err(|e| e.to_string())?;
            let current = state
                .rows
                .get(row)
                .and_then(|r| r.run_origin.as_ref())
                .ok_or("The saved origin changed.")?;
            if material_revision(current.as_ref())? != p.origin_revision {
                return Err("The saved origin changed.".into());
            }
            apply(
                state,
                row,
                r,
                &p.preview,
                &p.preview.context_revision,
                fragment,
            )
        }
        (InspectionRequest::Exact(r), InspectionPreview::Exact(p)) => {
            apply(state, row, r, p, revision, fragment)
        }
        (InspectionRequest::Mismatch(r), InspectionPreview::Mismatch(p)) => {
            let InspectionRequest::Mismatch(mut current) = build_inspection_request(state, row)?
            else {
                return Err("The PCR preparation mode changed.".into());
            };
            current.base.options = r.base.options.clone();
            let key = InspectionRequest::Mismatch(r.clone()).key()?;
            if InspectionRequest::Mismatch(current).key() != Ok(key.clone())
                || preview.request().key() != Ok(key)
            {
                return Err(
                    "The primers, mismatch policy, source role or catalog changed; inspect again."
                        .into(),
                );
            }
            p.validate_binding(revision).map_err(|e| e.to_string())?;
            let current_origin_revision = state.rows[row]
                .run_origin
                .as_ref()
                .map(|origin| material_revision(origin.as_ref()))
                .transpose()?;
            if current_origin_revision != p.origin_revision {
                return Err(
                    "The saved-material origin changed; inspect the mismatch PCR digest again."
                        .into(),
                );
            }
            if !preview.selectable() {
                return Err(
                    "The PCR geometry or physical digest fragments are unavailable.".into(),
                );
            }
            if !p
                .inspection
                .as_ref()
                .is_some_and(|i| i.fragments.iter().any(|f| f.index == fragment))
            {
                return Err("The selected fragment is unavailable.".into());
            }
            let number = fragment.checked_add(1).ok_or("Fragment number overflow.")?;
            let mut next = state.clone();
            next.rows[row].fragment = number.to_string();
            Ok(next)
        }
        (InspectionRequest::Extension(r), InspectionPreview::Extension(p)) => {
            let InspectionRequest::Extension(mut current) = build_inspection_request(state, row)?
            else {
                return Err("The overlap-extension preparation mode changed.".into());
            };
            current.options = r.options.clone();
            let key = InspectionRequest::Extension(r.clone()).key()?;
            if InspectionRequest::Extension(current).key() != Ok(key.clone())
                || preview.request().key() != Ok(key)
            {
                return Err("The six oligos, extension geometry, source role or catalog changed; inspect again.".into());
            }
            p.validate_binding(revision).map_err(|e| e.to_string())?;
            let current_origin_revision = state.rows[row]
                .run_origin
                .as_ref()
                .map(|origin| material_revision(origin.as_ref()))
                .transpose()?;
            if current_origin_revision != p.origin_revision {
                return Err("The saved-material origin changed; inspect the overlap-extension digest again.".into());
            }
            if !p.selectable() {
                return Err(
                    "The actual overlap-extension product or digest fragments are unavailable."
                        .into(),
                );
            }
            if !p.inspection.as_ref().is_some_and(|inspection| {
                inspection
                    .fragments
                    .iter()
                    .any(|item| item.index == fragment)
            }) {
                return Err("The selected outer-PCR digest fragment is unavailable.".into());
            }
            let number = fragment.checked_add(1).ok_or("Fragment number overflow.")?;
            let mut next = state.clone();
            next.rows[row].fragment = number.to_string();
            Ok(next)
        }
        _ => Err("The PCR inspection mode changed; inspect again.".into()),
    }
}

#[cfg(test)]
mod mismatch_digest_tests {
    use super::*;
    use collagene_core::cloning::{
        pcr::mismatch_digest as mismatch_pcr_digest_inspection,
        restriction::catalog as catalog_digest,
    };
    use collagene_model::{Molecule, Topology};
    const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTA";
    const MUTANT: &str = "AGACCGATGCTACGTAGGACTTCAGACCTGGA";
    fn fixture() -> (MethodsState, Molecule, InspectionRequest) {
        let mut state = MethodsState::default();
        let r = &mut state.rows[1];
        r.preparation = OverlapPreparation::MismatchPcrCatalogDigest;
        r.forward_name = "Mutant F + EcoRI".into();
        r.forward = "AGACCGAT".into();
        r.forward_tail = "CCGAATTC".into();
        r.reverse_name = "Mutant R + XhoI".into();
        r.reverse_primer = "TCCAGGTC".into();
        r.reverse_tail = "GGCTCGAG".into();
        r.mismatch_max_mismatches = "1".into();
        r.mismatch_three_prime_exact = "2".into();
        r.interval_start = "4".into();
        r.interval_end = "29".into();
        r.reverse = true;
        r.fragment = "19".into();
        let cat = catalog_digest::builtin_catalog().unwrap();
        let ids = cat
            .edition()
            .enzymes
            .iter()
            .filter(|e| e.name == "EcoRI" || e.name == "XhoI")
            .map(|e| e.id.clone())
            .collect::<Vec<_>>();
        r.catalog_selection = Some(catalog_digest::select(&cat, &ids).unwrap());
        let source = Molecule::dna("Original template", ORIGINAL, Topology::Linear).unwrap();
        let request = build_inspection_request(&state, 1).unwrap();
        (state, source, request)
    }
    fn inspect(source: &Molecule, request: &InspectionRequest) -> InspectionPreview {
        let InspectionRequest::Mismatch(r) = request else {
            panic!("expected mismatch request")
        };
        InspectionPreview::Mismatch(
            mismatch_pcr_digest_inspection::inspect(source, None, r, None).unwrap(),
        )
    }
    #[test]
    fn derived_mismatch_digest_keeps_native_origin_and_rejects_stale_application() {
        use collagene_core::{
            chemistry::Strand,
            documents::workspace::Workspace,
            material::{derivation as d, run::MaterialOrigin},
        };
        let (mut state, source, request) = fixture();
        let w = Workspace::from_documents("Retained mismatch template", vec![source.clone()], 0)
            .unwrap();
        let doc = &w.documents[0];
        let event = collagene_core::material::lineage::EventBinding {
            document_id: doc.id.clone(),
            event_id: doc.head_event.clone().unwrap(),
            revision: doc.head.clone(),
        };
        let forward = Primer::new("Template F", "AGTCCGAT", "").unwrap();
        let reverse = Primer::new("Template R", "TACAGGTC", "").unwrap();
        let original = d::SourceOrigin::OriginalReference;
        let documents =
            d::pcr::preview_pcr(&source, &original, &forward, &reverse, Strand::Top).unwrap();
        assert_eq!(documents.len(), 1);
        assert_eq!(documents[0].molecule.sequence.as_str(), ORIGINAL);
        let prepared = d::prepare_install(
            &w,
            &d::Request {
                sources: vec![d::BoundSource {
                    role: d::SourceRole::Template,
                    event: event,
                    origin: original.clone(),
                }],
                derivation: d::Derivation::ExactPcr {
                    selection: d::PcrSelection {
                        selected_index: 0,
                        expected_preview_revision: documents[0].molecule.revision(),
                        geometry: documents[0].geometry.clone(),
                    },
                    forward,
                    reverse,
                    template_strand: Strand::Top,
                },
                name: "Native mismatch-PCR template".into(),
                save_primers: false,
            },
        )
        .unwrap();
        let origin = MaterialOrigin::Derived(
            d::capsule(&prepared.workspace, &prepared.receipt.product).unwrap(),
        );
        let source = &prepared.workspace.revisions[&origin.event().revision];
        state.rows[1].run_origin = Some(
            crate::history::run_origin_model::PinnedOrigin::from_material(origin.clone()).unwrap(),
        );
        let InspectionRequest::Mismatch(native) = &request else {
            panic!("Mismatch")
        };
        let report = mismatch_pcr_digest_inspection::inspect(
            source,
            Some(&origin),
            native,
            Some(&prepared.workspace),
        )
        .unwrap();
        assert_eq!(
            report.amplicon.as_ref().unwrap().sequence.as_str(),
            format!("CCGAATTC{MUTANT}CTCGAGCC")
        );
        let revision = report.context_revision.clone();
        let preview = InspectionPreview::Mismatch(report);
        let next = apply_inspection(&state, 1, &request, &preview, &revision, 1).unwrap();
        assert!(next.rows[1].run_origin == state.rows[1].run_origin);
        assert_eq!(next.rows[1].fragment, "2");
        let mut changed = state.clone();
        changed.rows[1].forward_tail.push('A');
        assert!(apply_inspection(&changed, 1, &request, &preview, &revision, 1).is_err());
        state.rows[1].run_origin = None;
        assert!(apply_inspection(&state, 1, &request, &preview, &revision, 1).is_err());
    }
    #[test]
    fn actual_mutant_digest_selection_preserves_all_source_intent_and_oligo_settings() {
        let (state, source, request) = fixture();
        let preview = inspect(&source, &request);
        assert_eq!(
            preview.amplicon().unwrap().sequence.as_str(),
            format!("CCGAATTC{MUTANT}CTCGAGCC")
        );
        assert_eq!(preview.inspection().unwrap().fragments.len(), 3);
        assert!(!preview.incomplete());
        let revision = request.revision(&source, None, None).unwrap();
        let next = apply_inspection(&state, 1, &request, &preview, &revision, 1).unwrap();
        let mut expected = state.clone();
        expected.rows[1].fragment = "2".into();
        assert!(next == expected);
        let (_, rows) =
            crate::cloning::preparation::cloning_preparations::prepared_rows(&next, 2).unwrap();
        let collagene_core::cloning::assembly::recipe::Preparation::MismatchPcrCatalogDigest {
            options,
            forward,
            reverse,
            fragment_index,
            ..
        } = &rows[1].preparation
        else {
            panic!("selected preparation lost mismatch policy")
        };
        assert_eq!(options.binding.max_mismatches, 1);
        assert_eq!(options.binding.min_three_prime_exact, 2);
        assert_eq!(forward.ordered_sequence(), "CCGAATTCAGACCGAT");
        assert_eq!(reverse.ordered_sequence(), "GGCTCGAGTCCAGGTC");
        assert_eq!(*fragment_index, 1);
    }
    #[test]
    fn stale_policies_modes_and_fragments_cannot_reuse_a_mutant_receipt() {
        let (state, source, request) = fixture();
        let preview = inspect(&source, &request);
        let revision = request.revision(&source, None, None).unwrap();
        for field in 0..5 {
            let mut changed = state.clone();
            let row = &mut changed.rows[1];
            match field {
                0 => row.mismatch_max_mismatches = "0".into(),
                1 => row.mismatch_three_prime_exact = "3".into(),
                2 => row.preparation = OverlapPreparation::PcrCatalogDigest,
                3 => row.forward_name = "Different ordered oligo".into(),
                _ => row.reverse_tail.push('A'),
            }
            assert!(apply_inspection(&changed, 1, &request, &preview, &revision, 1).is_err());
        }
        assert!(apply_inspection(&state, 1, &request, &preview, &"f".repeat(64), 1).is_err());
        assert!(apply_inspection(&state, 1, &request, &preview, &revision, 99).is_err());
    }
    #[test]
    fn conditional_singleton_and_complete_amplicon_keep_the_recorded_mismatch_policy() {
        let (mut state, _, _) = fixture();
        state.rows[1].mismatch_max_products = "1".into();
        let request = build_inspection_request(&state, 1).unwrap();
        let repeated =
            Molecule::dna("Repeated original", &ORIGINAL.repeat(2), Topology::Linear).unwrap();
        let preview = inspect(&repeated, &request);
        assert!(preview.incomplete() && preview.selectable());
        let revision = request.revision(&repeated, None, None).unwrap();
        let applied = apply_inspection(&state, 1, &request, &preview, &revision, 1).unwrap();
        let (_, rows) =
            crate::cloning::preparation::cloning_preparations::prepared_rows(&applied, 2).unwrap();
        let collagene_core::cloning::assembly::recipe::Preparation::MismatchPcrCatalogDigest {
            options,
            ..
        } = &rows[1].preparation
        else {
            unreachable!()
        };
        assert_eq!(options.max_products, 1);
        let old = state.rows[1].clone();
        state.rows[1].use_complete_pcr_amplicon().unwrap();
        assert!(state.rows[1].preparation == OverlapPreparation::MismatchPcr);
        assert_eq!(
            state.rows[1].mismatch_max_mismatches,
            old.mismatch_max_mismatches
        );
        assert_eq!(
            state.rows[1].mismatch_three_prime_exact,
            old.mismatch_three_prime_exact
        );
        assert_eq!(state.rows[1].forward, old.forward);
        assert_eq!(state.rows[1].reverse_tail, old.reverse_tail);
    }
}

#[cfg(test)]
mod overlap_extension_digest_tests {
    use super::*;
    use crate::assembly::draft::assembly_role_ids::Role;
    use crate::cloning::preparation::cloning_methods_model::ExtensionPrimerInput;
    use crate::cloning::preparation::cloning_methods_model::Mode;
    use crate::cloning::preparation::cloning_methods_model::OverlapExtensionPcrInput;
    use crate::cloning::preparation::cloning_methods_model::Treatment;
    use collagene_core::{
        cloning::{overlap::extension::digest as native, restriction::catalog as catalog_digest},
        verification::evidence::Status,
    };
    use collagene_model::{Molecule, Topology};
    const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
    const MUTANT: &str = "AGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCA";
    fn oligo(name: &str, annealing: &str, tail: &str) -> ExtensionPrimerInput {
        ExtensionPrimerInput {
            name: name.into(),
            annealing: annealing.into(),
            tail: tail.into(),
        }
    }
    fn fixture() -> (MethodsState, Molecule, InspectionRequest) {
        let mut state = MethodsState::default();
        state.mode = Mode::Ligation;
        state.rows[0].preparation = OverlapPreparation::WholeLinear;
        state.source_roles = Some(vec![
            Role {
                role: "backbone".into(),
                input_id: "original-backbone".into(),
            },
            Role {
                role: "edited-insert".into(),
                input_id: "original-insert".into(),
            },
        ]);
        let row = &mut state.rows[1];
        row.preparation = OverlapPreparation::OverlapExtensionPcrCatalogDigest;
        row.document = 17;
        row.interval_start = "4".into();
        row.interval_end = "31".into();
        row.reverse = true;
        row.treatment = Treatment::Phosphorylate;
        row.repair = Some(collagene_core::chemistry::EndRepairMethod::Blunt);
        row.fragment = "19".into();
        row.right_fusion = "CGCG".into();
        row.forward = "Unused old single-pair field".into();
        let mut raw = OverlapExtensionPcrInput::default();
        raw.left.forward = oligo("Left outer", "AGTCCGAT", "");
        raw.left.reverse = oligo("Left mutagenic", "TACGTAGC", "GGTCTGAAAA");
        raw.right.forward = oligo("Right mutagenic", "TTCAGACC", "GCTACGTATT");
        raw.right.reverse = oligo("Right outer", "TGACGATC", "");
        raw.outer.forward = oligo("Outer F + EcoRI", "AGTCCGAT", "CCGAATTC");
        raw.outer.reverse = oligo("Outer R + XhoI", "TGACGATC", "GGCTCGAG");
        raw.overlap = "18".into();
        raw.min_paired_bases = "8".into();
        row.overlap_extension_pcr = Some(raw);
        let catalog = catalog_digest::builtin_catalog().unwrap();
        let ids = catalog
            .edition()
            .enzymes
            .iter()
            .filter(|entry| entry.name == "EcoRI" || entry.name == "XhoI")
            .map(|entry| entry.id.clone())
            .collect::<Vec<_>>();
        row.catalog_selection = Some(catalog_digest::select(&catalog, &ids).unwrap());
        let source =
            Molecule::dna("Unchanged original insert", ORIGINAL, Topology::Linear).unwrap();
        let request = build_inspection_request(&state, 1).unwrap();
        (state, source, request)
    }
    fn inspect(
        source: &Molecule,
        origin: Option<&collagene_core::material::run::MaterialOrigin>,
        request: &InspectionRequest,
        workspace: Option<&collagene_core::documents::workspace::Workspace>,
    ) -> InspectionPreview {
        let InspectionRequest::Extension(request) = request else {
            panic!("four-step inspection request");
        };
        InspectionPreview::Extension(native::inspect(source, origin, request, workspace).unwrap())
    }
    #[test]
    fn real_four_step_digest_selection_changes_only_fragment_and_retains_exact_source_and_oligos() {
        let (state, source, request) = fixture();
        let before = state.clone();
        let preview = inspect(&source, None, &request, None);
        let revision = request.revision(&source, None, None).unwrap();
        assert_eq!(
            preview.amplicon().unwrap().sequence.as_str(),
            format!("CCGAATTC{MUTANT}CTCGAGCC")
        );
        assert_eq!(preview.inspection().unwrap().fragments.len(), 3);
        assert_eq!(
            preview.inspection().unwrap().fragments[1].duplex.top_length,
            44
        );
        assert!(preview.selectable() && !preview.incomplete());
        assert!(preview.geometry_page(0).is_empty());
        let InspectionPreview::Extension(report) = &preview else {
            unreachable!()
        };
        assert_eq!(report.pcr_stages.len(), 3);
        assert_eq!(report.pcr_stages[0].template_revision, source.revision());
        assert_ne!(report.pcr_stages[2].template_revision, source.revision());
        let next = apply_inspection(&state, 1, &request, &preview, &revision, 1).unwrap();
        let mut expected = state.clone();
        expected.rows[1].fragment = "2".into();
        assert!(next == expected);
        assert!(state == before);
        let (_, rows) =
            crate::cloning::preparation::cloning_preparations::prepared_rows(&next, 2).unwrap();
        let collagene_core::cloning::assembly::recipe::Preparation::OverlapExtensionPcrCatalogDigest {
            conditions: _,
            request: prep,
            selection,
            fragment_index,
            repair,
        } = &rows[1].preparation
        else {
            panic!("preserved OE catalog chain");
        };
        assert_eq!(prep.left.reverse.ordered_sequence(), "GGTCTGAAAATACGTAGC");
        assert_eq!(prep.right.forward.ordered_sequence(), "GCTACGTATTTTCAGACC");
        assert_eq!(prep.outer.forward.tail, "CCGAATTC");
        assert_eq!(prep.outer.reverse.tail, "GGCTCGAG");
        assert_eq!(selection, next.rows[1].catalog_selection.as_ref().unwrap());
        assert_eq!(*fragment_index, 1);
        assert_eq!(*repair, next.rows[1].repair);
        assert_eq!(source.sequence.as_str(), ORIGINAL);
    }
    #[test]
    fn stale_any_of_six_oligos_or_geometry_catalog_role_and_missing_fragments_never_apply() {
        let (state, source, request) = fixture();
        let preview = inspect(&source, None, &request, None);
        let revision = request.revision(&source, None, None).unwrap();
        for field in 0..9 {
            let mut changed = state.clone();
            if field < 6 {
                let raw = changed.rows[1].overlap_extension_pcr.as_mut().unwrap();
                let primer = match field {
                    0 => &mut raw.left.forward,
                    1 => &mut raw.left.reverse,
                    2 => &mut raw.right.forward,
                    3 => &mut raw.right.reverse,
                    4 => &mut raw.outer.forward,
                    _ => &mut raw.outer.reverse,
                };
                primer.name.push_str(" changed");
            } else {
                match field {
                    6 => {
                        changed.rows[1]
                            .overlap_extension_pcr
                            .as_mut()
                            .unwrap()
                            .overlap = "17".into()
                    }
                    7 => changed.source_roles.as_mut().unwrap()[1]
                        .role
                        .push_str(" other"),
                    _ => changed.rows[1].catalog_selection = None,
                }
            }
            assert!(apply_inspection(&changed, 1, &request, &preview, &revision, 1).is_err());
        }
        assert!(apply_inspection(&state, 1, &request, &preview, &"0".repeat(64), 1).is_err());
        assert!(apply_inspection(&state, 1, &request, &preview, &revision, 99).is_err());
    }
    #[test]
    fn capped_extension_selection_preserves_the_cap_and_remains_incomplete_on_fresh_verification() {
        let (mut state, source, _) = fixture();
        state.rows[1]
            .overlap_extension_pcr
            .as_mut()
            .unwrap()
            .max_comparisons = "18".into();
        state.rows[1].repair = None;
        let request = build_inspection_request(&state, 1).unwrap();
        let preview = inspect(&source, None, &request, None);
        assert!(preview.incomplete() && preview.selectable());
        assert!(preview.amplicon().is_some());
        assert!(preview.inspection().is_some());
        let revision = request.revision(&source, None, None).unwrap();
        let next = apply_inspection(&state, 1, &request, &preview, &revision, 1).unwrap();
        assert_eq!(
            next.rows[1]
                .overlap_extension_pcr
                .as_ref()
                .unwrap()
                .max_comparisons,
            "18"
        );
        let (ids, fragments) =
            crate::cloning::preparation::cloning_preparations::prepared_rows(&next, 2).unwrap();
        let other =
            Molecule::dna("Other original source", "CCGGAATTGGCC", Topology::Linear).unwrap();
        let inputs = ids
            .into_iter()
            .zip([other, source])
            .map(|(id, molecule)| collagene_core::planning::plan::PlanInput {
                id,
                revision: molecule.revision(),
                molecule,
            })
            .collect::<Vec<_>>();
        let target = Molecule::dna(
            "Independent selected conditional digest",
            &format!("AATTC{MUTANT}C"),
            Topology::Linear,
        )
        .unwrap();
        let compiled = collagene_core::cloning::assembly::recipe::compile_preparations(
            &inputs,
            &target,
            "Fresh conditional preparation",
            &fragments,
            collagene_core::cloning::assembly::recipe::CleavagePolicy::AssumedComplete,
        )
        .unwrap();
        let plan = collagene_core::planning::plan::DesignPlan {
            name: "Fresh conditional preparation".into(),
            inputs,
            target_revision: target.revision(),
            steps: compiled.steps,
        };
        let checked = collagene_core::verification::engine::verify(&plan, &target);
        assert_eq!(checked.status, Status::Incomplete, "{:?}", checked.findings);
        assert_eq!(
            checked.products.last().unwrap().sequence.as_str(),
            target.sequence.as_str()
        );
    }
    #[test]
    fn saved_origin_is_replayed_and_cannot_be_removed_from_an_existing_digest_receipt() {
        use collagene_core::{
            documents::workspace::Workspace,
            material::{
                lineage::{EventBinding, InputEventBinding},
                outputs as cloning_outputs, run as run_material,
            },
            planning::plan::{DesignPlan, Operation, PlanInput, Source},
        };
        let (mut state, original, _) = fixture();
        let target = Molecule::dna(
            "Independent saved-stock PCR requirement",
            ORIGINAL,
            Topology::Linear,
        )
        .unwrap();
        let frame = Workspace::from_documents(
            "Actual original-stock preparation",
            vec![original.clone(), target.clone()],
            0,
        )
        .unwrap();
        let event = |index: usize| {
            let doc = &frame.documents[index];
            EventBinding {
                document_id: doc.id.clone(),
                event_id: doc.head_event.clone().unwrap(),
                revision: doc.head.clone(),
            }
        };
        let plan = DesignPlan {
            name: "Actual stock PCR".into(),
            target_revision: target.revision(),
            inputs: vec![PlanInput {
                id: "stock".into(),
                revision: original.revision(),
                molecule: original,
            }],
            steps: vec![Operation::Pcr {
                template_strand: collagene_core::chemistry::Strand::Top,
                template: Source::Input("stock".into()),
                forward: Primer::new("Stock F", "AGTCCGAT", "").unwrap(),
                reverse: Primer::new("Stock R", "TGACGATC", "").unwrap(),
            }],
        };
        let commit = cloning_outputs::finalize(
            &frame,
            &plan,
            event(1),
            vec![InputEventBinding {
                input_id: "stock".into(),
                event: event(0),
            }],
            &cloning_outputs::Options {
                name: "Saved original-stock PCR".into(),
                save_primers: false,
                show_prepared: false,
                allow_incomplete: false,
            },
        )
        .unwrap();
        let selection = run_material::Selection {
            run_id: commit.receipt.run.run_id.clone(),
            step_index: 0,
            event: commit.receipt.construct.clone(),
        };
        let origin =
            collagene_core::material::capsule::capsule(&commit.workspace, &selection).unwrap();
        let source = commit.workspace.revisions[&selection.event.revision].clone();
        state.rows[1].run_origin = Some(
            crate::history::run_origin_model::PinnedOrigin::from_material(
                collagene_core::material::run::MaterialOrigin::Run(origin.clone()),
            )
            .unwrap(),
        );
        let request = build_inspection_request(&state, 1).unwrap();
        assert!(request.revision(&source, None, None).is_err());
        let revision = request
            .revision(
                &source,
                Some(&collagene_core::material::run::MaterialOrigin::Run(
                    origin.clone(),
                )),
                Some(&commit.workspace),
            )
            .unwrap();
        let material = collagene_core::material::run::MaterialOrigin::Run(origin.clone());
        let preview = inspect(&source, Some(&material), &request, Some(&commit.workspace));
        let InspectionPreview::Extension(report) = &preview else {
            unreachable!()
        };
        assert_eq!(report.origin_incomplete, Some(false));
        assert_eq!(
            report.origin_revision.as_deref(),
            Some(
                material_revision(&collagene_core::material::run::MaterialOrigin::Run(
                    origin.clone()
                ))
                .unwrap()
                .as_str()
            )
        );
        let next = apply_inspection(&state, 1, &request, &preview, &revision, 1).unwrap();
        assert!(next.rows[1].run_origin == state.rows[1].run_origin);
        assert_eq!(next.rows[1].fragment, "2");
        state.rows[1].run_origin = None;
        assert!(apply_inspection(&state, 1, &request, &preview, &revision, 1).is_err());
    }
}
