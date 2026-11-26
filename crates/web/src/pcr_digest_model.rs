//! Pure PCR-product digest selection. Source events and UI epochs stay with the
//! caller; no original contribution interval is inferred from an amplicon.
use crate::cloning_methods_model::{MethodsState, OverlapPreparation};
use collagene_core::{
    analysis::Primer,
    pcr_digest_inspection::{Preview, Request},
    restriction_inspection::InspectionOptions,
};

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
    a.source_role == b.source_role
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
    use crate::{
        assembly_role_ids::Role,
        cloning_methods_model::{Mode, Treatment},
    };
    use collagene_core::{Molecule, Topology, catalog_digest, pcr_digest_inspection};

    fn select(name: &str) -> collagene_core::catalog_digest::CatalogSelection {
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

/// GUI-only sum types preserve both existing and additive native wire layouts.
#[derive(Clone, serde::Serialize)]
pub(crate) enum InspectionRequest {
    Exact(Request),
    Mismatch(collagene_core::mismatch_pcr_digest_inspection::Request),
}
pub(crate) enum InspectionPreview {
    Exact(Preview),
    Mismatch(collagene_core::mismatch_pcr_digest_inspection::Preview),
}
pub(crate) fn build_inspection_request(
    state: &MethodsState,
    row: usize,
) -> Result<InspectionRequest, String> {
    let base = build_base(state, row)?;
    if state.rows[row].preparation == OverlapPreparation::MismatchPcrCatalogDigest {
        Ok(InspectionRequest::Mismatch(
            collagene_core::mismatch_pcr_digest_inspection::Request {
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
        source: &collagene_core::Molecule,
        origin: Option<&collagene_core::run_material::Capsule>,
    ) -> Result<String, collagene_core::Error> {
        match self {
            Self::Exact(r) => match origin {
                Some(o) => collagene_core::pcr_digest_inspection::inspection_revision_from_run(
                    source, o, r,
                ),
                None => collagene_core::pcr_digest_inspection::inspection_revision(source, r),
            },
            Self::Mismatch(r) => {
                collagene_core::mismatch_pcr_digest_inspection::inspection_revision(
                    source, origin, r,
                )
            }
        }
    }
    pub(crate) fn into_job(
        self,
        molecule: collagene_core::Molecule,
        origin: Option<collagene_core::run_material::Capsule>,
    ) -> collagene_core::protocol::AnalysisJob {
        use collagene_core::protocol::AnalysisJob;
        match self {
            Self::Exact(request) => match origin {
                Some(origin) => AnalysisJob::InspectRunPcrDigest {
                    molecule,
                    origin,
                    request,
                },
                None => AnalysisJob::InspectPcrDigest { molecule, request },
            },
            Self::Mismatch(request) => AnalysisJob::InspectMismatchPcrDigest {
                molecule,
                origin,
                request,
            },
        }
    }
}
impl InspectionPreview {
    pub(crate) fn context_revision(&self) -> &str {
        match self {
            Self::Exact(p) => &p.context_revision,
            Self::Mismatch(p) => &p.context_revision,
        }
    }
    pub(crate) fn source_revision(&self) -> &str {
        match self {
            Self::Exact(p) => &p.source_revision,
            Self::Mismatch(p) => &p.source_revision,
        }
    }
    pub(crate) fn request(&self) -> InspectionRequest {
        match self {
            Self::Exact(p) => InspectionRequest::Exact(p.request.clone()),
            Self::Mismatch(p) => InspectionRequest::Mismatch(p.request.clone()),
        }
    }
    pub(crate) fn validate_binding(&self, r: &str) -> Result<(), collagene_core::Error> {
        match self {
            Self::Exact(p) => p.validate_binding(r),
            Self::Mismatch(p) => p.validate_binding(r),
        }
    }
    pub(crate) fn amplicon(&self) -> Option<&collagene_core::Molecule> {
        match self {
            Self::Exact(p) => p.amplicon.as_ref(),
            Self::Mismatch(p) => p.amplicon.as_ref(),
        }
    }
    pub(crate) fn inspection(
        &self,
    ) -> Option<&collagene_core::restriction_inspection::RestrictionInspection> {
        match self {
            Self::Exact(p) => p.inspection.as_ref(),
            Self::Mismatch(p) => p.inspection.as_ref(),
        }
    }
    pub(crate) fn findings(&self) -> &[String] {
        match self {
            Self::Exact(p) => &p.findings,
            Self::Mismatch(p) => &p.findings,
        }
    }
    pub(crate) fn geometry_count(&self) -> usize {
        match self {
            Self::Exact(p) => p.pcr_geometries.len(),
            Self::Mismatch(p) => p.pcr.report.products.len(),
        }
    }
    pub(crate) fn geometry_page(
        &self,
        start: usize,
    ) -> Vec<collagene_core::pcr_digest_inspection::GeometrySummary> {
        match self {
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
                    |(index, g)| collagene_core::pcr_digest_inspection::GeometrySummary {
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
        matches!(self,Self::Mismatch(p) if p.incomplete)
    }
    pub(crate) fn selectable(&self) -> bool {
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
        _ => Err("The PCR inspection mode changed; inspect again.".into()),
    }
}

#[cfg(test)]
mod mismatch_digest_tests {
    use super::*;
    use collagene_core::{Molecule, Topology, catalog_digest, mismatch_pcr_digest_inspection};
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
            mismatch_pcr_digest_inspection::inspect(source, None, r).unwrap(),
        )
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
        let revision = request.revision(&source, None).unwrap();
        let next = apply_inspection(&state, 1, &request, &preview, &revision, 1).unwrap();
        let mut expected = state.clone();
        expected.rows[1].fragment = "2".into();
        assert!(next == expected);
        let (_, rows) = crate::cloning_preparations::prepared_rows(&next, 2).unwrap();
        let collagene_core::assembly_recipe::Preparation::MismatchPcrCatalogDigest {
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
        let revision = request.revision(&source, None).unwrap();
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
        let revision = request.revision(&repeated, None).unwrap();
        let applied = apply_inspection(&state, 1, &request, &preview, &revision, 1).unwrap();
        let (_, rows) = crate::cloning_preparations::prepared_rows(&applied, 2).unwrap();
        let collagene_core::assembly_recipe::Preparation::MismatchPcrCatalogDigest {
            options, ..
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
