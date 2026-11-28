use super::model::*;
use crate::{
    chemistry::tracked_duplex::{self, CheckedDigestOptions, CleavageBudget, TrackedDuplex},
    cloning::{pcr::source::Template, restriction::catalog},
    encoding::bounded_json,
    material::{
        derivation::{self, SourceOrigin},
        run::MaterialOrigin,
    },
    sequence::analysis::ExactSearchBudget,
    verification::replay_budget::ReplayBudget,
};
use collagene_model::{Error, Topology};
use std::collections::BTreeSet;

fn bad(message: &str) -> Error {
    Error::Format(message.into())
}
fn text(value: &str, maximum: usize) -> bool {
    !value.is_empty() && value.len() <= maximum && !value.chars().any(char::is_control)
}
fn concrete(value: &str) -> bool {
    value.bytes().all(|b| b"ACGT".contains(&b))
}
impl Request {
    pub fn validate(&self) -> Result<(), Error> {
        self.limits.validate()?;
        if self.components.is_empty()
            || self.components.len() > MAX_COMPONENTS
            || !(2..=MAX_SUBJECTS).contains(&self.subjects.len())
            || self.tests.is_empty()
            || self.tests.len() > MAX_TESTS
            || self.readout.counterpart_tolerance_twice_bp > 2_000_000
        {
            return Err(Error::TooLarge);
        }
        // Admit all opaque packages and model trees before any closed-origin
        // decode, hashing or origin replay. Share the canonical mixed ledger.
        let mut opaque_bytes = 0usize;
        for component in &self.components {
            collagene_model::bounds::model_preflight(&[&component.molecule])?;
            let bytes = match &component.origin {
                SourceOrigin::OriginalReference => 0,
                SourceOrigin::Saved(MaterialOrigin::Run(origin)) => origin.candidate.len(),
                SourceOrigin::Saved(MaterialOrigin::Derived(origin)) => origin.record_bytes.len(),
            };
            opaque_bytes = opaque_bytes
                .checked_add(bytes)
                .filter(|bytes| *bytes <= MAX_REQUEST_BYTES)
                .ok_or(Error::TooLarge)?;
        }
        crate::material::capsule::material_origins_preflight(self.components.iter().filter_map(
            |component| {
                if let SourceOrigin::Saved(origin) = &component.origin {
                    Some(origin)
                } else {
                    None
                }
            },
        ))?;
        bounded_json::preflight(self, MAX_REQUEST_BYTES)?;
        let mut ids = BTreeSet::new();
        let mut bases = 0usize;
        for component in &self.components {
            component.molecule.validate()?;
            bases = bases
                .checked_add(component.molecule.sequence.len())
                .ok_or(Error::TooLarge)?;
            if bases > 1_000_000 {
                return Err(Error::TooLarge);
            }
            if !text(&component.id, 96)
                || !text(&component.label, 256)
                || !ids.insert(component.id.as_str())
                || component.revision != component.molecule.revision()
            {
                return Err(bad("invalid component identity or exact source revision"));
            }
            if matches!(&component.origin, SourceOrigin::Saved(_)) {
                if component.original_generation.is_some() {
                    return Err(bad(
                        "saved material cannot be relabeled as original generation",
                    ));
                }
            }
            if let Some(generation) = &component.original_generation {
                generation.validate()?;
            }
        }
        let mut subject_ids = BTreeSet::new();
        let mut intended = 0;
        for subject in &self.subjects {
            if !text(&subject.id, 96)
                || !text(&subject.label, 256)
                || !subject_ids.insert(subject.id.as_str())
                || subject.components.is_empty()
                || subject.components.len() > MAX_COMPONENTS
            {
                return Err(bad("invalid named subject or component list"));
            }
            intended += usize::from(subject.role == Role::Intended);
            let mut components = BTreeSet::new();
            for id in &subject.components {
                if !ids.contains(id.as_str()) || !components.insert(id.as_str()) {
                    return Err(bad("unknown or duplicated subject component"));
                }
            }
        }
        if intended != 1 {
            return Err(bad(
                "exactly one intended subject and at least one named decoy required",
            ));
        }
        // Unused component snapshots are refused: every retained source has an explicit subject.
        if self
            .components
            .iter()
            .any(|c| !self.subjects.iter().any(|s| s.components.contains(&c.id)))
        {
            return Err(bad("component is not associated with any named subject"));
        }
        let current = catalog::builtin_catalog()?;
        let mut test_ids = BTreeSet::new();
        for test in &self.tests {
            if !text(&test.id, 96) || !text(&test.label, 256) || !test_ids.insert(test.id.as_str())
            {
                return Err(bad("invalid or duplicate screening test identity"));
            }
            match &test.method {
                Method::Digest { selection } => {
                    selection.validate_for(&current)?;
                    if selection.enzymes.len() > 4 {
                        return Err(Error::TooLarge);
                    }
                }
                Method::ExactPcr {
                    forward, reverse, ..
                } => {
                    for primer in [forward, reverse] {
                        if !text(&primer.name, 256)
                            || primer.annealing.is_empty()
                            || primer.annealing.len() > 200
                            || primer.tail.len() > 200
                            || !primer.tail.bytes().all(|b| b"ACGTRYSWKMBDHVN".contains(&b))
                        {
                            return Err(bad("invalid bounded PCR primer"));
                        }
                    }
                }
            }
        }
        Ok(())
    }
}
pub fn revision(request: &Request) -> Result<String, Error> {
    request.validate()?;
    bounded_json::hash(&(ENGINE, request))
}

struct Material {
    tracked: TrackedDuplex,
    revision: String,
    incomplete: bool,
}
fn material(component: &Component, replay: &mut ReplayBudget<'_>) -> Result<Material, Error> {
    // Canonical original guards or one shared native saved-origin replay run
    // before any override. Nested exact-PCR/cleavage cannot reset this request.
    let (mut tracked, incomplete) = match &component.origin {
        SourceOrigin::OriginalReference => {
            let fresh = derivation::fresh_source(&component.molecule, &component.origin)?;
            let (_, tracked, incomplete, _, _, _) = fresh.into_parts();
            (tracked, incomplete)
        }
        SourceOrigin::Saved(origin) => {
            let (fresh, _) = origin.replay_scoped(&component.molecule, replay, None)?;
            let (_, tracked, incomplete, _, _, _) = fresh.into_parts();
            (tracked, incomplete)
        }
    };
    if let Some(generation) = &component.original_generation {
        tracked = TrackedDuplex::from_input(&component.molecule, generation)?;
    }
    let revision = tracked.summary()?.state_revision;
    Ok(Material {
        tracked,
        revision,
        incomplete,
    })
}
fn blank(test_index: usize, component_index: usize, component: &Component) -> Cell {
    Cell {
        test_index,
        component_index,
        source_revision: component.revision.clone(),
        material_revision: None,
        state: CellState::Unsearched,
        sizes_twice_bp: vec![],
        native_product_count: 0,
        enumeration_complete: false,
        pcr_geometries: vec![],
        cleavage: None,
        reason: String::new(),
    }
}
fn cell(
    request: &Request,
    test_index: usize,
    component_index: usize,
    material: &Result<Material, Error>,
    cleavage: &mut CleavageBudget,
    pcr: &mut ExactSearchBudget,
) -> Result<Cell, Error> {
    let component = &request.components[component_index];
    let mut result = blank(test_index, component_index, component);
    let material = match material {
        Ok(material) => material,
        Err(error) => {
            result.state = if matches!(error, Error::TooLarge)
                || pcr.remaining() == 0
                || cleavage.remaining() == 0
            {
                CellState::Limited
            } else {
                CellState::Unsupported
            };
            result.reason = error.to_string().chars().take(1024).collect();
            return Ok(result);
        }
    };
    result.material_revision = Some(material.revision.clone());
    if material.incomplete {
        result.state = CellState::Ambiguous;
        result.reason =
            "Fresh source replay has unresolved evidence; no screening absence claim.".into();
        return Ok(result);
    }
    match &request.tests[test_index].method {
        Method::Digest { selection } => {
            if cleavage.remaining() == 0 {
                result.reason = "Shared native cleavage budget exhausted before this cell.".into();
                return Ok(result);
            }
            let digest = tracked_duplex::checked_digest_with_budget(
                &material.tracked,
                &CheckedDigestOptions {
                    selection: selection.clone(),
                    max_comparisons: request.limits.max_cleavage_work,
                    max_sites: request.limits.max_sites_per_cell,
                },
                cleavage,
            );
            let digest = match digest {
                Ok(value) => value,
                Err(Error::TooLarge) => {
                    result.state = CellState::Limited;
                    result.reason =
                        "Native checked-cleavage resource bound refused this cell.".into();
                    return Ok(result);
                }
                Err(error) => {
                    result.state = CellState::Unsupported;
                    result.reason = error.to_string().chars().take(1024).collect();
                    return Ok(result);
                }
            };
            result.enumeration_complete = digest.complete && digest.products_determined;
            result.native_product_count = digest.fragments.len();
            result.cleavage = Some(digest.evidence());
            if !digest.complete {
                result.state = CellState::Limited;
                result.reason =
                    "Incomplete native recognition/context scan; missing cuts not excluded.".into();
            } else if !digest.products_determined {
                result.state = CellState::Ambiguous;
                result.reason =
                    "Positional cleavage remains unresolved; no determined fragment signature."
                        .into();
            } else if digest.fragments.len() > request.limits.max_products_per_cell {
                result.state = CellState::Limited;
                result.reason = "Native fragment count exceeds retained-signature bound.".into();
            } else if digest
                .fragments
                .iter()
                .any(|f| f.duplex.duplex.topology != Topology::Linear)
            {
                result.state = CellState::Unsupported;
                result.reason = "Uncut circular material has no supported linear-fragment size signature; it is not no signal.".into();
            } else {
                result.state = CellState::Complete;
                result.sizes_twice_bp = digest
                    .fragments
                    .iter()
                    .map(|f| f.duplex.duplex.top.len() + f.duplex.duplex.bottom.len())
                    .collect();
            }
        }
        Method::ExactPcr {
            forward,
            reverse,
            template_strand,
        } => {
            if pcr.remaining() == 0 {
                result.reason = "Shared native exact-PCR budget exhausted before this cell.".into();
                return Ok(result);
            }
            let template =
                Template::fresh(&component.molecule, &material.tracked, *template_strand);
            let products = template.and_then(|t| t.geometries_with_budget(forward, reverse, pcr));
            match products {
                Ok((products, _materialized)) => {
                    result.enumeration_complete = true;
                    result.native_product_count = products.len();
                    if products.len() > request.limits.max_products_per_cell {
                        result.state = CellState::Limited;
                        result.reason = "Complete finite native geometry enumeration exceeds retained-product bound; no absence/distinction claim.".into();
                    } else {
                        let d = &material.tracked.duplex;
                        let certain = concrete(d.top.as_str())
                            && concrete(d.bottom.as_str())
                            && concrete(forward.annealing.as_str())
                            && concrete(reverse.annealing.as_str())
                            && concrete(&forward.tail)
                            && concrete(&reverse.tail);
                        result.state = if certain {
                            CellState::Complete
                        } else {
                            CellState::Ambiguous
                        };
                        result.reason = if certain {
                            "Finite exact native PCR geometry enumeration; products are modeled possibilities, not experimental amplification claims.".into()
                        } else {
                            "IUPAC ambiguity prevents a definitive size/absence signature.".into()
                        };
                        for geometry in products {
                            let product = geometry.product;
                            result.sizes_twice_bp.push(2 * product.sequence.len());
                            result.pcr_geometries.push(PcrGeometry {
                                forward_start: product.forward_start,
                                reverse_start: product.reverse_start,
                                crosses_origin: product.crosses_origin,
                                product_length: product.sequence.len(),
                            });
                        }
                    }
                }
                Err(error) => {
                    result.state = if pcr.remaining() == 0 || matches!(error, Error::TooLarge) {
                        CellState::Limited
                    } else {
                        CellState::Unsupported
                    };
                    result.reason = error.to_string().chars().take(1024).collect();
                }
            }
        }
    }
    result.sizes_twice_bp.sort_unstable();
    result.sizes_twice_bp.dedup();
    Ok(result)
}

/// A size must lack ANY nearby counterpart. Different counts or greedy
/// one-to-one matching do not establish separation without intensity assumptions.
pub(super) fn distinguish(intended: &[usize], decoy: &[usize], tolerance: usize) -> Distinction {
    // Both native signatures are sorted distinct sets. Advance the opposite
    // pointer monotonically instead of rescanning up to 1024 sizes for each size.
    for (sizes, opposite, side) in [
        (intended, decoy, WitnessSide::Intended),
        (decoy, intended, WitnessSide::Decoy),
    ] {
        let mut counterpart = 0usize;
        for &size in sizes {
            let lower = size.saturating_sub(tolerance);
            while counterpart < opposite.len() && opposite[counterpart] < lower {
                counterpart += 1;
            }
            if opposite
                .get(counterpart)
                .is_none_or(|other| size.abs_diff(*other) > tolerance)
            {
                return Distinction::Separated(Witness {
                    side,
                    size_twice_bp: size,
                });
            }
        }
    }
    if intended == decoy {
        Distinction::IdenticalSignature
    } else {
        Distinction::WithinDeclaredTolerance
    }
}
fn merge_state(a: &CellState, b: &CellState) -> CellState {
    // Every component must be complete before a mixture can certify absence.
    for state in [
        CellState::Unsearched,
        CellState::Limited,
        CellState::Unsupported,
        CellState::Ambiguous,
    ] {
        if a == &state || b == &state {
            return state;
        }
    }
    CellState::Complete
}
fn test_result(
    subjects: &[Subject],
    component_ids: &[&str],
    definition: &Test,
    test_index: usize,
    cells: &[Cell],
    readout: &IdealSizeReadout,
) -> TestResult {
    let signatures: Vec<_> = subjects
        .iter()
        .enumerate()
        .map(|(subject_index, subject)| {
            let mut state = CellState::Complete;
            let mut sizes = BTreeSet::new();
            for id in &subject.components {
                let index = component_ids
                    .iter()
                    .position(|c| *c == id.as_str())
                    .expect("validated component");
                let cell = &cells[index];
                state = merge_state(&state, &cell.state);
                sizes.extend(cell.sizes_twice_bp.iter().copied());
            }
            Signature {
                subject_index,
                state,
                sizes_twice_bp: sizes.into_iter().collect(),
            }
        })
        .collect();
    let intended = subjects
        .iter()
        .position(|s| s.role == Role::Intended)
        .expect("validated intended");
    let comparisons: Vec<_> = subjects
        .iter()
        .enumerate()
        .filter(|(_, s)| s.role != Role::Intended)
        .map(|(index, _)| {
            let a = &signatures[intended];
            let b = &signatures[index];
            let distinction = if a.state == CellState::Complete && b.state == CellState::Complete {
                distinguish(
                    &a.sizes_twice_bp,
                    &b.sizes_twice_bp,
                    readout.counterpart_tolerance_twice_bp,
                )
            } else {
                Distinction::Unresolved
            };
            Comparison {
                decoy_subject_index: index,
                distinction,
            }
        })
        .collect();
    TestResult {
        test_index,
        definition: definition.clone(),
        separated_decoys: comparisons
            .iter()
            .filter(|c| matches!(c.distinction, Distinction::Separated(_)))
            .count(),
        complete_comparisons: comparisons
            .iter()
            .filter(|c| c.distinction != Distinction::Unresolved)
            .count(),
        signatures,
        comparisons,
    }
}

pub(super) fn panel_search(tests: &[TestResult], maximum: usize) -> PanelSearch {
    let mut result = PanelSearch {
        test_indices: vec![],
        separated_decoy_indices: vec![],
        evaluated_combinations: 0,
        comparison_work: 0,
        complete_through_size_three: true,
    };
    // 64 choose 1+2+3 <= 43744. Work charged for each visited comparison,
    // including repeated/shared decoys, not just retained candidate panels.
    let mut work = 0;
    'search: for size in 1..=3 {
        for a in 0..tests.len() {
            for b in a..tests.len() {
                if size == 1 && b != a {
                    break;
                }
                if size >= 2 && b == a {
                    continue;
                }
                for c in b..tests.len() {
                    if size <= 2 && c != b {
                        break;
                    }
                    if size == 3 && c == b {
                        continue;
                    }
                    let indices = match size {
                        1 => vec![a],
                        2 => vec![a, b],
                        _ => vec![a, b, c],
                    };
                    let charge = 1 + indices
                        .iter()
                        .map(|&i| tests[i].comparisons.len())
                        .sum::<usize>();
                    if charge > maximum.saturating_sub(work) {
                        result.complete_through_size_three = false;
                        break 'search;
                    }
                    work += charge;
                    let mut separated = BTreeSet::new();
                    for &index in &indices {
                        for comparison in &tests[index].comparisons {
                            if matches!(comparison.distinction, Distinction::Separated(_)) {
                                separated.insert(comparison.decoy_subject_index);
                            }
                        }
                    }
                    result.evaluated_combinations += 1;
                    if separated.len() > result.separated_decoy_indices.len() {
                        result.test_indices = indices;
                        result.separated_decoy_indices = separated.into_iter().collect();
                    }
                }
            }
        }
    }
    result.comparison_work = work;
    result
}

fn json_bytes(value: &impl serde::Serialize, limit: usize) -> Result<usize, Error> {
    struct Count {
        used: usize,
        limit: usize,
    }
    impl std::io::Write for Count {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.used = self
                .used
                .checked_add(bytes.len())
                .filter(|n| *n <= self.limit)
                .ok_or_else(|| std::io::Error::other("screening report retention bound"))?;
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut count = Count { used: 0, limit };
    serde_json::to_writer(&mut count, value).map_err(|_| Error::TooLarge)?;
    Ok(count.used)
}

pub fn assess(request: &Request) -> Result<Report, Error> {
    let context_revision = revision(request)?;
    bounded_json::preflight(
        &(
            &request.subjects,
            &request.tests,
            &request.readout,
            &request.limits,
            SCOPE,
        ),
        request.limits.max_report_bytes,
    )?;
    let mut cleavage = CleavageBudget::new(request.limits.max_cleavage_work)?;
    let mut replay = ReplayBudget::shared_job(&mut cleavage);
    replay.pcr = ExactSearchBudget::new(request.limits.max_pcr_work)?;
    let materials: Vec<_> = request
        .components
        .iter()
        .map(|component| material(component, &mut replay))
        .collect();
    let mut report = Report {
        engine_version: ENGINE.into(),
        context_revision,
        readout: request.readout.clone(),
        limits: request.limits.clone(),
        components: request
            .components
            .iter()
            .map(|c| {
                Ok(ComponentBinding {
                    id: c.id.clone(),
                    label: c.label.clone(),
                    source_revision: c.revision.clone(),
                    origin: match &c.origin {
                        SourceOrigin::OriginalReference => ComponentOrigin::OriginalReference,
                        SourceOrigin::Saved(_) => ComponentOrigin::SavedNative,
                    },
                    origin_revision: bounded_json::hash(&(
                        "screening-source-1",
                        &c.origin,
                        &c.original_generation,
                    ))?,
                    original_generation: c.original_generation.clone(),
                })
            })
            .collect::<Result<Vec<_>, Error>>()?,
        subjects: request.subjects.clone(),
        cells: vec![],
        tests: vec![],
        ranked_test_indices: vec![],
        panel: PanelSearch {
            test_indices: vec![],
            separated_decoy_indices: vec![],
            evaluated_combinations: 0,
            comparison_work: 0,
            complete_through_size_three: false,
        },
        cleavage_work: 0,
        pcr_work: 0,
        statement: SCOPE.into(),
        report_revision: String::new(),
    };
    let mut retained_cell_bytes = 0usize;
    for test_index in 0..request.tests.len() {
        let start = report.cells.len();
        for component_index in 0..request.components.len() {
            let cell = cell(
                request,
                test_index,
                component_index,
                &materials[component_index],
                replay.cleavage,
                &mut replay.pcr,
            )?;
            // Bound retained cell and the aggregate before publishing it to the report.
            retained_cell_bytes = retained_cell_bytes
                .checked_add(json_bytes(
                    &cell,
                    request.limits.max_report_bytes - retained_cell_bytes,
                )?)
                .ok_or(Error::TooLarge)?;
            report.cells.push(cell);
        }
        report.tests.push(test_result(
            &request.subjects,
            &request
                .components
                .iter()
                .map(|c| c.id.as_str())
                .collect::<Vec<_>>(),
            &request.tests[test_index],
            test_index,
            &report.cells[start..],
            &request.readout,
        ));
    }
    report.ranked_test_indices = (0..report.tests.len()).collect();
    report.ranked_test_indices.sort_by_key(|&index| {
        let result = &report.tests[index];
        (
            std::cmp::Reverse(result.separated_decoys),
            std::cmp::Reverse(result.complete_comparisons),
            index,
        )
    });
    report.panel = panel_search(&report.tests, request.limits.max_panel_checks);
    report.cleavage_work = replay.cleavage_used();
    report.pcr_work = replay.pcr_used();
    // Reserve the final identity's fixed-size text before hashing.
    report.report_revision = "0".repeat(64);
    bounded_json::preflight(&report, request.limits.max_report_bytes)?;
    report.report_revision = report.identity()?;
    Ok(report)
}
impl Report {
    fn identity(&self) -> Result<String, Error> {
        bounded_json::hash(&(
            &self.engine_version,
            &self.context_revision,
            &self.readout,
            &self.limits,
            &self.components,
            &self.subjects,
            &self.cells,
            &self.tests,
            &self.ranked_test_indices,
            &self.panel,
            self.cleavage_work,
            self.pcr_work,
            &self.statement,
        ))
    }
    /// Documentary identity/shape only, never a substitute for fresh native replay.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.limits.validate()?;
        bounded_json::preflight(self, self.limits.max_report_bytes)?;
        if self.engine_version != ENGINE
            || self.context_revision != expected
            || self.statement != SCOPE
            || self.report_revision != self.identity()?
            || !(2..=MAX_SUBJECTS).contains(&self.subjects.len())
            || self.tests.is_empty()
            || self.tests.len() > MAX_TESTS
            || self.components.is_empty()
            || self.components.len() > MAX_COMPONENTS
            || self.cells.len() != self.components.len() * self.tests.len()
            || self.ranked_test_indices.len() != self.tests.len()
            || self.panel.test_indices.len() > 3
            || self.panel.comparison_work > self.limits.max_panel_checks
            || self.cleavage_work > self.limits.max_cleavage_work
            || self.pcr_work > self.limits.max_pcr_work
        {
            return Err(bad("screening report binding/shape mismatch"));
        }
        let mut ranked = self.ranked_test_indices.clone();
        ranked.sort_unstable();
        if ranked != (0..self.tests.len()).collect::<Vec<_>>() {
            return Err(bad(
                "screening ranking must retain each original test index",
            ));
        }
        if self.context_revision.len() != 64
            || self.report_revision.len() != 64
            || self.readout.counterpart_tolerance_twice_bp > 2_000_000
            || self
                .subjects
                .iter()
                .filter(|s| s.role == Role::Intended)
                .count()
                != 1
        {
            return Err(bad("invalid screening documentary context"));
        }
        let mut component_ids = BTreeSet::new();
        for component in &self.components {
            if !text(&component.id, 96)
                || !text(&component.label, 256)
                || !component_ids.insert(component.id.as_str())
                || component.source_revision.len() != 64
                || component.origin_revision.len() != 64
                || (component.origin == ComponentOrigin::SavedNative
                    && component.original_generation.is_some())
            {
                return Err(bad("invalid screening source descriptor"));
            }
            if let Some(generation) = &component.original_generation {
                generation.validate()?;
            }
        }
        let mut subject_ids = BTreeSet::new();
        for subject in &self.subjects {
            if !text(&subject.id, 96)
                || !text(&subject.label, 256)
                || !subject_ids.insert(subject.id.as_str())
                || subject.components.is_empty()
                || subject.components.len() > MAX_COMPONENTS
            {
                return Err(bad("invalid documentary named subject"));
            }
            let mut ids = BTreeSet::new();
            for id in &subject.components {
                if !component_ids.contains(id.as_str()) || !ids.insert(id.as_str()) {
                    return Err(bad("invalid documentary mixture association"));
                }
            }
        }
        for (index, test) in self.tests.iter().enumerate() {
            if test.test_index != index
                || test.signatures.len() != self.subjects.len()
                || test.comparisons.len() != self.subjects.len() - 1
            {
                return Err(bad("screening test/subject slots differ"));
            }
        }
        for (slot, cell) in self.cells.iter().enumerate() {
            if cell.test_index != slot / self.components.len()
                || cell.component_index != slot % self.components.len()
                || cell.source_revision != self.components[cell.component_index].source_revision
                || cell.sizes_twice_bp.len() > 64
                || cell.pcr_geometries.len() > 64
                || cell.reason.len() > 4096
                || cell.source_revision.len() != 64
                || cell
                    .material_revision
                    .as_ref()
                    .is_some_and(|r| r.len() != 64)
                || (cell.state == CellState::Complete && !cell.enumeration_complete)
                || cell
                    .sizes_twice_bp
                    .windows(2)
                    .any(|pair| pair[0] >= pair[1])
            {
                return Err(bad("invalid screening cell"));
            }
            if let Some(evidence) = &cell.cleavage {
                evidence.validate()?;
                if cell.material_revision.as_deref() != Some(evidence.source_revision.as_str())
                    || (cell.enumeration_complete
                        && (!evidence.complete || !evidence.products_determined))
                {
                    return Err(bad(
                        "cleavage evidence differs from actual material/enumeration binding",
                    ));
                }
            }
            if cell.state == CellState::Complete
                && (cell.native_product_count > self.limits.max_products_per_cell
                    || cell.material_revision.is_none())
            {
                return Err(bad("complete screening cell lacks bounded native evidence"));
            }
            match &self.tests[cell.test_index].definition.method {
                Method::ExactPcr { .. } => {
                    if cell.cleavage.is_some()
                        || (cell.state == CellState::Complete
                            && cell.native_product_count != cell.pcr_geometries.len())
                    {
                        return Err(bad("PCR screening cell has inconsistent evidence"));
                    }
                    let mut sizes = Vec::with_capacity(cell.pcr_geometries.len());
                    for geometry in &cell.pcr_geometries {
                        if geometry.product_length == 0
                            || geometry.product_length > collagene_model::MAX_SEQUENCE_LEN
                            || geometry.forward_start >= collagene_model::MAX_SEQUENCE_LEN
                            || geometry.reverse_start >= collagene_model::MAX_SEQUENCE_LEN
                        {
                            return Err(bad(
                                "PCR screening geometry is outside canonical reference bounds",
                            ));
                        }
                        sizes.push(
                            geometry
                                .product_length
                                .checked_mul(2)
                                .ok_or(Error::TooLarge)?,
                        );
                    }
                    sizes.sort_unstable();
                    sizes.dedup();
                    if cell.sizes_twice_bp != sizes {
                        return Err(bad(
                            "PCR size signature differs from retained native geometry lengths",
                        ));
                    }
                }
                Method::Digest { .. } => {
                    if !cell.pcr_geometries.is_empty()
                        || (cell.state == CellState::Complete
                            && (cell.cleavage.is_none()
                                || cell.native_product_count == 0
                                || cell.sizes_twice_bp.is_empty()))
                    {
                        return Err(bad("digest screening cell has inconsistent evidence"));
                    }
                }
            }
        }
        let ids: Vec<_> = self.components.iter().map(|c| c.id.as_str()).collect();
        for (index, test) in self.tests.iter().enumerate() {
            let start = index * self.components.len();
            let expected = test_result(
                &self.subjects,
                &ids,
                &test.definition,
                index,
                &self.cells[start..start + self.components.len()],
                &self.readout,
            );
            if test.signatures != expected.signatures
                || test.comparisons != expected.comparisons
                || test.separated_decoys != expected.separated_decoys
                || test.complete_comparisons != expected.complete_comparisons
            {
                return Err(bad(
                    "screening conclusions differ from retained component signatures",
                ));
            }
        }
        let mut expected_ranking: Vec<_> = (0..self.tests.len()).collect();
        expected_ranking.sort_by_key(|&i| {
            (
                std::cmp::Reverse(self.tests[i].separated_decoys),
                std::cmp::Reverse(self.tests[i].complete_comparisons),
                i,
            )
        });
        if self.ranked_test_indices != expected_ranking
            || self.panel != panel_search(&self.tests, self.limits.max_panel_checks)
        {
            return Err(bad(
                "screening rank/panel differs from bounded declared comparison policy",
            ));
        }
        Ok(())
    }
    pub fn validate_for(&self, request: &Request) -> Result<(), Error> {
        self.validate_binding(&revision(request)?)?;
        if self.report_revision != assess(request)?.report_revision {
            return Err(bad("screening report differs from fresh native assessment"));
        }
        Ok(())
    }
}
