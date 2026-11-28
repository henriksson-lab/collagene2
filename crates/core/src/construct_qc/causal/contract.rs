//! Bounded current request admission and documentary report bindings.
use super::*;
use crate::{
    encoding::bounded_json,
    material::{lineage, run},
    sequence::translation,
};
use std::collections::BTreeSet;

fn invalid(message: &str) -> Error {
    Error::Format(format!("causal construct QC: {message}"))
}
fn pin(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}
fn text(value: &str, max: usize) -> bool {
    !value.trim().is_empty() && value.len() <= max && !value.chars().any(char::is_control)
}
fn assumption_shape(assumption: &crate::construct_qc::AnnotationAssumption) -> Result<(), Error> {
    if assumption.feature_index >= 4096
        || !pin(&assumption.source_revision)
        || !pin(&assumption.annotation_revision)
        || !text(&assumption.rationale, 1024)
    {
        return Err(invalid(
            "bounded pinned annotation-role assumptions required",
        ));
    }
    Ok(())
}
fn transition_shape(cause: &TransitionPin) -> Result<(), Error> {
    match cause {
        TransitionPin::ReferenceEdit {
            event,
            source,
            actions_revision,
        } => {
            run::event_shape(event)?;
            run::event_shape(source)?;
            if !pin(actions_revision) {
                return Err(invalid("invalid native edit context"));
            }
        }
        TransitionPin::Fork { event, source } => {
            run::event_shape(event)?;
            run::event_shape(source)?;
            if event.document_id == source.document_id || event.revision != source.revision {
                return Err(invalid("invalid fork identity receipt"));
            }
        }
        TransitionPin::CloningStep {
            event,
            run_id,
            step_index,
            plan_revision,
        } => {
            run::event_shape(event)?;
            if !pin(run_id) || !pin(plan_revision) || *step_index >= 1000 {
                return Err(invalid("invalid native run product context"));
            }
        }
        TransitionPin::MaterialDerivation {
            event,
            context_revision,
            record_hash,
        } => {
            run::event_shape(event)?;
            if !pin(context_revision) || !pin(record_hash) {
                return Err(invalid("invalid native material context"));
            }
        }
        TransitionPin::TopoGraph {
            record_hash,
            state_revision,
            context_revision,
        } => {
            if !pin(record_hash) || !pin(state_revision) || !pin(context_revision) {
                return Err(invalid("invalid native graph context"));
            }
        }
    }
    Ok(())
}
fn codon_location(location: &collagene_model::Location) -> Result<(), Error> {
    use collagene_model::annotations::biological_path::{BiologicalPath, PathLimits};
    let path = BiologicalPath::from_location(
        location,
        collagene_model::MAX_SEQUENCE_LEN,
        PathLimits {
            max_calls: 3,
            max_spans: 3,
        },
    )?;
    if path.len() != 3 {
        return Err(invalid(
            "codon evidence requires exactly three resolved biological calls",
        ));
    }
    Ok(())
}
impl Request {
    pub fn validate(&self, workspace: &Workspace) -> Result<(), Error> {
        self.validate_shape()?;
        workspace.validate()
    }
    pub(super) fn validate_shape(&self) -> Result<(), Error> {
        if self.schema_version != 1 {
            return Err(invalid("unsupported request schema"));
        }
        let limits = &self.limits;
        if self.pairs.is_empty()
            || self.pairs.len() > 32
            || !(1..=1000).contains(&limits.max_history_events)
            || !(1..=2_000_000).contains(&limits.max_coding_calls)
            || !(1..=50_000_000).contains(&limits.max_work)
            || !(1..=512).contains(&limits.max_segments)
            || !(1..=256).contains(&limits.max_stops)
        {
            return Err(Error::TooLarge);
        }
        let mut ids = BTreeSet::new();
        for pair in &self.pairs {
            if !text(&pair.id, 128) || !ids.insert(&pair.id) {
                return Err(invalid("distinct bounded pair identities required"));
            }
            run::event_shape(&pair.source)?;
            run::event_shape(&pair.product)?;
            translation::CodeDescriptor::for_id(pair.genetic_code)?;
            for assumption in [&pair.source_cds, &pair.product_cds] {
                assumption_shape(assumption)?;
            }
        }
        bounded_json::preflight(self, MAX_REQUEST_BYTES)
    }
}
pub(super) fn revision(workspace: &Workspace, request: &Request) -> Result<String, Error> {
    request.validate_shape()?;
    let fingerprint = lineage::workspace_fingerprint(workspace)?;
    context(&fingerprint, request)
}
pub(super) fn request_revision(request: &Request) -> Result<String, Error> {
    request.validate_shape()?;
    bounded_json::hash(&(ENGINE_VERSION, request))
}
pub(super) fn context(fingerprint: &str, request: &Request) -> Result<String, Error> {
    bounded_json::hash(&(ENGINE_VERSION, fingerprint, request_revision(request)?))
}
pub(super) fn digest(report: &Report) -> Result<String, Error> {
    bounded_json::hash(&(
        &report.engine_version,
        &report.context_revision,
        &report.workspace_fingerprint,
        &report.request_revision,
        &report.pairs,
        report.work,
        report.coding_calls,
        report.history_events,
        &report.statement,
    ))
}
impl Report {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.validate_shape()?;
        if self.context_revision != expected || self.report_revision != digest(self)? {
            return Err(invalid(
                "report content or current native request binding changed",
            ));
        }
        Ok(())
    }
    /// Bounded current request/documentary validation, without scientific replay.
    pub fn validate_request(&self, workspace: &Workspace, request: &Request) -> Result<(), Error> {
        self.validate_binding(&revision(workspace, request)?)?;
        if self.request_revision != request_revision(request)?
            || self.pairs.len() != request.pairs.len()
            || self.work > request.limits.max_work
            || self.coding_calls > request.limits.max_coding_calls
            || self.history_events > request.limits.max_history_events
        {
            return Err(invalid("report does not belong to exact bounded request"));
        }
        let mut segments = 0usize;
        let mut stops = 0usize;
        for (actual, expected) in self.pairs.iter().zip(&request.pairs) {
            if actual.id != expected.id
                || actual.source != expected.source
                || actual.product != expected.product
                || actual.source_cds != expected.source_cds
                || actual.product_cds != expected.product_cds
                || actual.code.descriptor.id != expected.genetic_code
            {
                return Err(invalid("report pair declarations changed"));
            }
            if let (Ok(source), Ok(product)) = (
                run::resolve(workspace, &actual.source),
                run::resolve(workspace, &actual.product),
            ) {
                for finding in &actual.findings {
                    for location in &finding.product_coordinates {
                        location.validate(product.sequence.len())?;
                    }
                }
                for stop in &actual.new_stops {
                    stop.product_coordinates.validate(product.sequence.len())?;
                    if let Some(location) = &stop.parent_coordinates {
                        location.validate(source.sequence.len())?;
                    }
                }
            } else if actual.mapping_complete
                || !actual.register_segments.is_empty()
                || !actual.new_stops.is_empty()
            {
                return Err(invalid(
                    "unavailable event cannot carry complete coordinate evidence",
                ));
            }
            segments += actual.register_segments.len();
            stops += actual.new_stops.len();
        }
        if segments > request.limits.max_segments || stops > request.limits.max_stops {
            return Err(Error::TooLarge);
        }
        Ok(())
    }
    pub fn validate_for(&self, workspace: &Workspace, request: &Request) -> Result<(), Error> {
        self.validate_request(workspace, request)?;
        let actual = assess(workspace, request)?;
        if self != &actual {
            return Err(invalid(
                "report differs from fresh current native causal assessment",
            ));
        }
        Ok(())
    }
    pub(super) fn validate_shape(&self) -> Result<(), Error> {
        if self.engine_version != ENGINE_VERSION
            || self.statement != STATEMENT
            || !pin(&self.context_revision)
            || !pin(&self.workspace_fingerprint)
            || !pin(&self.request_revision)
            || !pin(&self.report_revision)
            || self.pairs.is_empty()
            || self.pairs.len() > 32
            || self.work > 50_000_000
            || self.coding_calls > 2_000_000
            || self.history_events > 1000
        {
            return Err(invalid("unsupported or oversized documentary report"));
        }
        let mut ids = BTreeSet::new();
        let mut locations = collagene_model::bounds::ModelBudget::default();
        let mut segments = 0usize;
        let mut stops = 0usize;
        for pair in &self.pairs {
            if !text(&pair.id, 128)
                || !ids.insert(&pair.id)
                || pair.findings.len() > 32
                || !text(&pair.annotation_correspondence.reason, 4096)
            {
                return Err(Error::TooLarge);
            }
            run::event_shape(&pair.source)?;
            run::event_shape(&pair.product)?;
            assumption_shape(&pair.source_cds)?;
            assumption_shape(&pair.product_cds)?;
            if pair.code != crate::construct_qc::CodePin::for_id(pair.code.descriptor.id)? {
                return Err(invalid("genetic code pin changed"));
            }
            segments = segments
                .checked_add(pair.register_segments.len())
                .ok_or(Error::TooLarge)?;
            stops = stops
                .checked_add(pair.new_stops.len())
                .ok_or(Error::TooLarge)?;
            if segments > 512 || stops > 256 {
                return Err(Error::TooLarge);
            }
            for finding in &pair.findings {
                if !text(&finding.reason, 4096)
                    || finding.product_coordinates.len() > 16
                    || finding.citations.is_empty()
                    || finding.citations.len() > 8
                {
                    return Err(invalid("invalid named causal finding"));
                }
                if matches!(finding.rule, Rule::JunctionRegister | Rule::NewInternalStop)
                    && matches!(
                        finding.verdict,
                        crate::construct_qc::Verdict::Pass | crate::construct_qc::Verdict::Fail
                    )
                    && !finding.under_declared_annotation_assumptions
                {
                    return Err(invalid(
                        "causal coding conclusions require explicit current role assumptions",
                    ));
                }
                for citation in &finding.citations {
                    if !text(&citation.source, 512) || !text(&citation.locator, 512) {
                        return Err(Error::TooLarge);
                    }
                }
                for location in &finding.product_coordinates {
                    locations.location(location)?;
                }
            }
            for segment in &pair.register_segments {
                if segment.product_coding_start >= segment.product_coding_end
                    || segment.product_coding_end > 2_000_000
                    || segment.first_parent_coding_ordinal >= 2_000_000
                    || segment.last_parent_coding_ordinal >= 2_000_000
                    || segment.offset_mod3 > 2
                    || segment.causes.len() > 64
                    || segment.native_junction_proven && segment.causes.is_empty()
                {
                    return Err(invalid("invalid register evidence"));
                }
                for cause in &segment.causes {
                    transition_shape(cause)?;
                }
            }
            for stop in &pair.new_stops {
                let has_parent = stop.parent_codon_index.is_some();
                if stop.product_codon_index >= 666_667
                    || stop.product_codon.len() != 3
                    || !stop
                        .product_codon
                        .bytes()
                        .all(|base| crate::sequence::analysis::base_mask(base) != 0)
                    || stop
                        .parent_codon_index
                        .is_some_and(|index| index >= 666_667)
                    || stop.parent_codon.as_ref().is_some_and(|codon| {
                        codon.len() != 3
                            || !codon
                                .bytes()
                                .all(|base| crate::sequence::analysis::base_mask(base) != 0)
                    })
                    || has_parent != stop.parent_coordinates.is_some()
                    || has_parent != stop.parent_codon.is_some()
                    || has_parent != stop.parent_amino_acid.is_some()
                    || stop
                        .parent_amino_acid
                        .is_some_and(|amino| !amino.is_ascii_uppercase() && amino != '*')
                    || !text(&stop.reason, 4096)
                    || stop.causes.len() > 64
                {
                    return Err(invalid("invalid stop evidence"));
                }
                for cause in &stop.causes {
                    transition_shape(cause)?;
                }
                locations.location(&stop.product_coordinates)?;
                codon_location(&stop.product_coordinates)?;
                if let Some(location) = &stop.parent_coordinates {
                    locations.location(location)?;
                    codon_location(location)?;
                }
            }
        }
        bounded_json::preflight(self, MAX_REPORT_BYTES)
    }
}
