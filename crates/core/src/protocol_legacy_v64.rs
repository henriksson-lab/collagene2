//! Exact pre-66 positional report layouts. Historical evidence gains no processed
//! graph and no current scientific authority during conversion.
use crate::{Error, Molecule, protocol, verification};
use serde::Deserialize;
const HEADER: usize = 42;
#[derive(Deserialize)]
struct Report {
    plan_revision: String,
    target_revision: String,
    engine_version: String,
    status: verification::Status,
    findings: Vec<String>,
    assumptions: Vec<String>,
    products: Vec<Molecule>,
    exhaustive_alternative_search: bool,
    step_reports: Vec<verification::StepReport>,
    annotation_reports: Vec<verification::AnnotationStepReport>,
}
impl From<Report> for verification::VerificationReport {
    fn from(r: Report) -> Self {
        Self {
            plan_revision: r.plan_revision,
            target_revision: r.target_revision,
            engine_version: r.engine_version,
            status: r.status,
            findings: r.findings,
            assumptions: r.assumptions,
            products: r.products,
            exhaustive_alternative_search: r.exhaustive_alternative_search,
            step_reports: r.step_reports,
            annotation_reports: r.annotation_reports,
            processed_step_reports: Vec::new(),
            physical_processing_reports: Vec::new(),
            reference_step_reports: Vec::new(),
        }
    }
}
#[derive(Deserialize)]
struct ReportV65 {
    plan_revision: String,
    target_revision: String,
    engine_version: String,
    status: verification::Status,
    findings: Vec<String>,
    assumptions: Vec<String>,
    products: Vec<Molecule>,
    exhaustive_alternative_search: bool,
    step_reports: Vec<verification::StepReport>,
    annotation_reports: Vec<verification::AnnotationStepReport>,
    processed_step_reports: Vec<verification::ProcessedTopoStepReport>,
}
impl From<ReportV65> for verification::VerificationReport {
    fn from(r: ReportV65) -> Self {
        Self {
            plan_revision: r.plan_revision,
            target_revision: r.target_revision,
            engine_version: r.engine_version,
            status: r.status,
            findings: r.findings,
            assumptions: r.assumptions,
            products: r.products,
            exhaustive_alternative_search: r.exhaustive_alternative_search,
            step_reports: r.step_reports,
            annotation_reports: r.annotation_reports,
            processed_step_reports: r.processed_step_reports,
            physical_processing_reports: Vec::new(),
            reference_step_reports: Vec::new(),
        }
    }
}
#[derive(Deserialize)]
struct ReportV67 {
    plan_revision: String,
    target_revision: String,
    engine_version: String,
    status: verification::Status,
    findings: Vec<String>,
    assumptions: Vec<String>,
    products: Vec<Molecule>,
    exhaustive_alternative_search: bool,
    step_reports: Vec<verification::StepReport>,
    annotation_reports: Vec<verification::AnnotationStepReport>,
    processed_step_reports: Vec<verification::ProcessedTopoStepReport>,
    physical_processing_reports: Vec<verification::PhysicalProcessingStepReport>,
}
impl From<ReportV67> for verification::VerificationReport {
    fn from(r: ReportV67) -> Self {
        Self {
            plan_revision: r.plan_revision,
            target_revision: r.target_revision,
            engine_version: r.engine_version,
            status: r.status,
            findings: r.findings,
            assumptions: r.assumptions,
            products: r.products,
            exhaustive_alternative_search: r.exhaustive_alternative_search,
            step_reports: r.step_reports,
            annotation_reports: r.annotation_reports,
            processed_step_reports: r.processed_step_reports,
            physical_processing_reports: r.physical_processing_reports,
            reference_step_reports: Vec::new(),
        }
    }
}
// Nested inverse-PCR and replay existed before physical/annotation reports.
// Their outer positional layouts stayed unchanged; retain each real report era.
#[derive(Deserialize)]
struct ReportV24 {
    plan_revision: String,
    target_revision: String,
    engine_version: String,
    status: verification::Status,
    findings: Vec<String>,
    assumptions: Vec<String>,
    products: Vec<Molecule>,
    exhaustive_alternative_search: bool,
}
impl From<ReportV24> for verification::VerificationReport {
    fn from(r: ReportV24) -> Self {
        Self {
            plan_revision: r.plan_revision,
            target_revision: r.target_revision,
            engine_version: r.engine_version,
            status: r.status,
            findings: r.findings,
            assumptions: r.assumptions,
            products: r.products,
            exhaustive_alternative_search: r.exhaustive_alternative_search,
            step_reports: Vec::new(),
            annotation_reports: Vec::new(),
            processed_step_reports: Vec::new(),
            physical_processing_reports: Vec::new(),
            reference_step_reports: Vec::new(),
        }
    }
}
#[derive(Deserialize)]
struct ReportV26 {
    plan_revision: String,
    target_revision: String,
    engine_version: String,
    status: verification::Status,
    findings: Vec<String>,
    assumptions: Vec<String>,
    products: Vec<Molecule>,
    exhaustive_alternative_search: bool,
    step_reports: Vec<verification::StepReport>,
}
impl From<ReportV26> for verification::VerificationReport {
    fn from(r: ReportV26) -> Self {
        Self {
            plan_revision: r.plan_revision,
            target_revision: r.target_revision,
            engine_version: r.engine_version,
            status: r.status,
            findings: r.findings,
            assumptions: r.assumptions,
            products: r.products,
            exhaustive_alternative_search: r.exhaustive_alternative_search,
            step_reports: r.step_reports,
            annotation_reports: Vec::new(),
            processed_step_reports: Vec::new(),
            physical_processing_reports: Vec::new(),
            reference_step_reports: Vec::new(),
        }
    }
}
// Each old outer struct retains fields *following* its nested report. Generic
// EOF/default migration cannot recover these boundaries in postcard.
#[derive(Deserialize)]
struct Inverse<R = Report> {
    source_revision: String,
    target_revision: String,
    options: crate::mutagenesis::InversePcrOptions,
    forward: crate::analysis::Primer,
    reverse: crate::analysis::Primer,
    forward_location: crate::Location,
    reverse_location: crate::Location,
    retained_template_length: usize,
    plan: verification::DesignPlan,
    verification: R,
    warnings: Vec<String>,
}
impl<R: Into<verification::VerificationReport>> From<Inverse<R>>
    for crate::mutagenesis::InversePcrDesign
{
    fn from(r: Inverse<R>) -> Self {
        Self {
            source_revision: r.source_revision,
            target_revision: r.target_revision,
            options: r.options,
            forward: r.forward,
            reverse: r.reverse,
            forward_location: r.forward_location,
            reverse_location: r.reverse_location,
            retained_template_length: r.retained_template_length,
            plan: r.plan,
            verification: r.verification.into(),
            warnings: r.warnings,
        }
    }
}
#[derive(Deserialize)]
struct Replay<R = Report> {
    run_id: String,
    source_workspace_fingerprint: String,
    recorded_status: verification::Status,
    recorded_engine: String,
    current_report: R,
    engine_matches: bool,
    products_match: bool,
    status_matches: bool,
}
impl<R: Into<verification::VerificationReport>> From<Replay<R>>
    for crate::run_lineage::ReplayReport
{
    fn from(r: Replay<R>) -> Self {
        Self {
            run_id: r.run_id,
            source_workspace_fingerprint: r.source_workspace_fingerprint,
            recorded_status: r.recorded_status,
            recorded_engine: r.recorded_engine,
            current_report: r.current_report.into(),
            engine_matches: r.engine_matches,
            products_match: r.products_match,
            status_matches: r.status_matches,
        }
    }
}
#[derive(Deserialize)]
struct Overlap<R = Report> {
    context_revision: String,
    compiled: crate::assembly_recipe::Compiled,
    verification: R,
    prepared: Vec<Molecule>,
    profile: Option<crate::overlap_profiles::ProfileReport>,
    preparation_complete: bool,
}
impl<R: Into<verification::VerificationReport>> From<Overlap<R>>
    for crate::assembly_recipe::Preview
{
    fn from(r: Overlap<R>) -> Self {
        Self {
            context_revision: r.context_revision,
            compiled: r.compiled,
            verification: r.verification.into(),
            prepared: r.prepared,
            profile: r.profile,
            preparation_complete: r.preparation_complete,
        }
    }
}
#[derive(Deserialize)]
struct Ligation<R = Report> {
    context_revision: String,
    compiled: crate::ligation_recipe::Compiled,
    verification: R,
    prepared: Vec<Molecule>,
    physical: Vec<Option<crate::chemistry::DuplexSummary>>,
    product_physical: Option<crate::chemistry::DuplexSummary>,
    junctions: Vec<crate::ligation_recipe::Junction>,
    preparation_complete: bool,
}
impl<R: Into<verification::VerificationReport>> From<Ligation<R>>
    for crate::ligation_recipe::Preview
{
    fn from(r: Ligation<R>) -> Self {
        Self {
            context_revision: r.context_revision,
            compiled: r.compiled,
            verification: r.verification.into(),
            prepared: r.prepared,
            physical: r.physical,
            product_physical: r.product_physical,
            junctions: r.junctions,
            preparation_complete: r.preparation_complete,
        }
    }
}
#[derive(Deserialize)]
struct Junction<R = Report> {
    engine_version: String,
    context_revision: String,
    plan_revision: String,
    target_revision: String,
    options: crate::assembly_junctions::Options,
    verification: R,
    product_revision: Option<String>,
    product_length: Option<usize>,
    product_topology: Option<crate::Topology>,
    method: Option<crate::assembly_junctions::Method>,
    junctions: Vec<crate::assembly_junctions::Junction>,
    complete: bool,
    omitted_junctions: usize,
    warnings: Vec<String>,
}
impl<R: Into<verification::VerificationReport>> From<Junction<R>>
    for crate::assembly_junctions::Report
{
    fn from(r: Junction<R>) -> Self {
        Self {
            engine_version: r.engine_version,
            context_revision: r.context_revision,
            plan_revision: r.plan_revision,
            target_revision: r.target_revision,
            options: r.options,
            verification: r.verification.into(),
            product_revision: r.product_revision,
            product_length: r.product_length,
            product_topology: r.product_topology,
            method: r.method,
            junctions: r.junctions,
            complete: r.complete,
            omitted_junctions: r.omitted_junctions,
            warnings: r.warnings,
        }
    }
}
#[derive(Deserialize)]
struct Pool<R = Report> {
    context_revision: String,
    request: crate::prepared_golden_gate::Request,
    target: Molecule,
    preview: Ligation<R>,
    missing_physical: Vec<usize>,
    pool: Option<crate::golden_gate::GoldenGateReport>,
    assessments: Vec<crate::prepared_golden_gate::Assessment>,
    warnings: Vec<String>,
}
impl<R: Into<verification::VerificationReport>> From<Pool<R>>
    for crate::prepared_golden_gate::Report
{
    fn from(r: Pool<R>) -> Self {
        Self {
            context_revision: r.context_revision,
            request: r.request,
            target: r.target,
            preview: r.preview.into(),
            missing_physical: r.missing_physical,
            pool: r.pool,
            assessments: r.assessments,
            warnings: r.warnings,
        }
    }
}
#[derive(Deserialize)]
struct Reaction<R = Report> {
    context_revision: String,
    request: crate::golden_gate_reaction::Request,
    target: Molecule,
    probe_plan: verification::DesignPlan,
    materials: Vec<crate::golden_gate_reaction::Material>,
    verification: R,
    pool: Option<crate::golden_gate::GoldenGateReport>,
    assessments: Vec<crate::prepared_golden_gate::Assessment>,
    warnings: Vec<String>,
}
impl<R: Into<verification::VerificationReport>> From<Reaction<R>>
    for crate::golden_gate_reaction::Report
{
    fn from(r: Reaction<R>) -> Self {
        Self {
            context_revision: r.context_revision,
            request: r.request,
            target: r.target,
            probe_plan: r.probe_plan,
            materials: r.materials,
            verification: r.verification.into(),
            pool: r.pool,
            assessments: r.assessments,
            warnings: r.warnings,
        }
    }
}
#[derive(Deserialize)]
struct Domestication<R = Report> {
    context_revision: String,
    request: crate::domestication_preparation::Request,
    source_revision: String,
    target: Molecule,
    check: crate::domestication::CheckReport,
    preparation: crate::mutagenesis::InversePcrPreparation,
    verification: R,
    warnings: Vec<String>,
}
impl<R: Into<verification::VerificationReport>> From<Domestication<R>>
    for crate::domestication_preparation::Prepared
{
    fn from(r: Domestication<R>) -> Self {
        Self {
            context_revision: r.context_revision,
            request: r.request,
            source_revision: r.source_revision,
            target: r.target,
            check: r.check,
            preparation: r.preparation,
            verification: r.verification.into(),
            warnings: r.warnings,
        }
    }
}
fn invalid(s: &str) -> Error {
    Error::Format(format!("historical report envelope: {s}"))
}
fn payload<T: for<'de> Deserialize<'de>>(bytes: &[u8]) -> Result<T, Error> {
    let (value, rest) = postcard::take_from_bytes(bytes).map_err(|e| invalid(&e.to_string()))?;
    if !rest.is_empty() {
        return Err(invalid("unexpected trailing payload"));
    }
    Ok(value)
}
/// The caller retains older1..26 bare-report converters and applies result
/// minimum versions and readonly response guards after this layout conversion.
pub(crate) fn decode(bytes: &[u8], maximum: usize) -> Result<Option<protocol::Response>, Error> {
    if bytes.len() < HEADER || &bytes[..4] != b"CGR\0" {
        return Err(invalid("unsupported magic/header"));
    }
    let version = u16::from_le_bytes([bytes[4], bytes[5]]);
    if version == 0 || version > protocol::WORKER_SCHEMA_VERSION {
        return Err(Error::Unsupported(format!(
            "binary schema version {version}"
        )));
    }
    let length = u32::from_le_bytes(bytes[6..10].try_into().unwrap()) as usize;
    if length > maximum {
        return Err(Error::TooLarge);
    }
    if HEADER.checked_add(length) != Some(bytes.len()) {
        return Err(invalid("binary length mismatch"));
    }
    if blake3::hash(&bytes[HEADER..]).as_bytes() != &bytes[10..42] {
        return Err(invalid("binary checksum mismatch"));
    }
    if version > 67 {
        return Ok(None);
    }
    let ((id, revision, tag), rest) =
        postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[HEADER..])
            .map_err(|e| invalid(&e.to_string()))?;
    if revision.len() > 128 {
        return Err(invalid("revision bound"));
    }
    use protocol::AnalysisResult as A;
    let result = match tag {
        4 if version >= 66 => A::Verification(payload::<ReportV67>(rest)?.into()),
        16 if version >= 66 => A::InversePcr(payload::<Inverse<ReportV67>>(rest)?.into()),
        41 if version >= 66 => A::RunReplay(payload::<Replay<ReportV67>>(rest)?.into()),
        86 if version >= 66 => A::OverlapRecipePreview(payload::<Overlap<ReportV67>>(rest)?.into()),
        87 if version >= 66 => {
            A::LigationRecipePreview(payload::<Ligation<ReportV67>>(rest)?.into())
        }
        88 if version >= 66 => {
            A::AssemblyJunctionInspection(payload::<Junction<ReportV67>>(rest)?.into())
        }
        105 if version >= 66 => A::PreparedGoldenGatePool(payload::<Pool<ReportV67>>(rest)?.into()),
        107 if version >= 66 => A::GoldenGateReaction(payload::<Reaction<ReportV67>>(rest)?.into()),
        134 if version >= 66 => {
            A::DomesticationPrepared(payload::<Domestication<ReportV67>>(rest)?.into())
        }
        4 if version == 65 => A::Verification(payload::<ReportV65>(rest)?.into()),
        16 if version == 65 => A::InversePcr(payload::<Inverse<ReportV65>>(rest)?.into()),
        41 if version == 65 => A::RunReplay(payload::<Replay<ReportV65>>(rest)?.into()),
        86 if version == 65 => A::OverlapRecipePreview(payload::<Overlap<ReportV65>>(rest)?.into()),
        87 if version == 65 => {
            A::LigationRecipePreview(payload::<Ligation<ReportV65>>(rest)?.into())
        }
        88 if version == 65 => {
            A::AssemblyJunctionInspection(payload::<Junction<ReportV65>>(rest)?.into())
        }
        105 if version == 65 => A::PreparedGoldenGatePool(payload::<Pool<ReportV65>>(rest)?.into()),
        107 if version == 65 => A::GoldenGateReaction(payload::<Reaction<ReportV65>>(rest)?.into()),
        134 if version == 65 => {
            A::DomesticationPrepared(payload::<Domestication<ReportV65>>(rest)?.into())
        }
        4 if version >= 27 => A::Verification(payload::<Report>(rest)?.into()),
        16 if version < 25 => A::InversePcr(payload::<Inverse<ReportV24>>(rest)?.into()),
        16 if version < 27 => A::InversePcr(payload::<Inverse<ReportV26>>(rest)?.into()),
        16 => A::InversePcr(payload::<Inverse>(rest)?.into()),
        41 if version < 25 => A::RunReplay(payload::<Replay<ReportV24>>(rest)?.into()),
        41 if version < 27 => A::RunReplay(payload::<Replay<ReportV26>>(rest)?.into()),
        41 => A::RunReplay(payload::<Replay>(rest)?.into()),
        86 => A::OverlapRecipePreview(payload::<Overlap>(rest)?.into()),
        87 => A::LigationRecipePreview(payload::<Ligation>(rest)?.into()),
        88 => A::AssemblyJunctionInspection(payload::<Junction>(rest)?.into()),
        105 => A::PreparedGoldenGatePool(payload::<Pool>(rest)?.into()),
        107 => A::GoldenGateReaction(payload::<Reaction>(rest)?.into()),
        134 => A::DomesticationPrepared(payload::<Domestication>(rest)?.into()),
        _ => return Ok(None),
    };
    Ok(Some(protocol::Response {
        id,
        revision: revision.into(),
        result,
    }))
}
