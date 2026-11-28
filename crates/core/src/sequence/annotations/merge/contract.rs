//! Borrowed resource admission and documentary report pins; no source replay here.
use super::*;
use crate::documents::workspace::Workspace;
use crate::material::{lineage::EventBinding, run};
use std::collections::BTreeSet;
fn library_preflight(
    library: &library::FeatureLibrary,
    report: &library::ScanReport,
) -> Result<(), Error> {
    library.validate()?;
    report.options.validate()?;
    if report.hits.len() > library::MAX_SCAN_HITS {
        return Err(Error::TooLarge);
    }
    let mut budget = collagene_model::bounds::ModelBudget::annotation_edit();
    for hit in &report.hits {
        budget.location(&hit.location)?;
        budget.features(std::slice::from_ref(&hit.feature))?;
        if !text(&hit.id, 128) {
            return Err(bad("Invalid bounded library proposal identifier"));
        }
    }
    bounded_json::preflight(&(library, report), MAX_REQUEST_BYTES - 1024)
}
fn dataset_preflight(
    report: &crate::formats::annotations::AnnotationImportReport,
) -> Result<(), Error> {
    if report.proposals.len() > crate::formats::annotations::MAX_ANNOTATION_ROWS
        || report.original_base64.len()
            > crate::formats::annotations::MAX_ANNOTATION_BYTES.div_ceil(3) * 4
    {
        return Err(Error::TooLarge);
    }
    let mut budget = collagene_model::bounds::ModelBudget::annotation_edit();
    for proposal in &report.proposals {
        budget.features(std::slice::from_ref(&proposal.feature))?;
        if !text(&proposal.id, 128) {
            return Err(bad("Invalid bounded dataset proposal identifier"));
        }
    }
    bounded_json::preflight(report, MAX_REQUEST_BYTES - 1024)
}
impl Source {
    pub fn from_library(
        library: &library::FeatureLibrary,
        report: &library::ScanReport,
    ) -> Result<Self, Error> {
        library_preflight(library, report)?;
        Ok(Self::Library {
            library: library.clone(),
            report: report.clone(),
        })
    }
    pub fn from_dataset(
        report: &crate::formats::annotations::AnnotationImportReport,
    ) -> Result<Self, Error> {
        dataset_preflight(report)?;
        Ok(Self::Dataset {
            report: report.clone(),
        })
    }

    /// Structural proposal identifiers only, not permission to apply source claims.
    pub fn proposal_ids(&self) -> Result<Vec<String>, Error> {
        self.preflight()?;
        Ok(match self {
            Self::Library { report, .. } => report.hits.iter().map(|p| p.id.clone()).collect(),
            Self::Dataset { report } => report.proposals.iter().map(|p| p.id.clone()).collect(),
            Self::Donor { indices, .. } => indices
                .iter()
                .map(|index| format!("donor-{index}"))
                .collect(),
        })
    }
    pub fn preflight(&self) -> Result<(), Error> {
        match self {
            Self::Library { library, report } => library_preflight(library, report)?,
            Self::Dataset { report } => dataset_preflight(report)?,
            Self::Donor {
                molecule,
                event,
                indices,
                ..
            } => {
                collagene_model::bounds::native_frame_preflight(
                    std::iter::once(molecule),
                    MAX_REQUEST_BYTES,
                )?;
                molecule.validate()?;
                if let Self::Donor { transform, .. } = self {
                    if molecule.sequence.alphabet() != collagene_model::Alphabet::Dna
                        || transform.origin >= molecule.sequence.len()
                        || (molecule.topology == collagene_model::Topology::Linear
                            && transform.origin != 0)
                    {
                        return Err(bad(
                            "Direct donor requires DNA and a valid declared origin transform",
                        ));
                    }
                }
                run::event_shape(event)?;
                if molecule.revision() != event.revision {
                    return Err(bad("Donor snapshot differs from declared event content"));
                }
                if indices.is_empty()
                    || indices.len() > MAX_SELECTIONS
                    || indices.iter().any(|i| *i >= molecule.features.len())
                    || indices.iter().collect::<BTreeSet<_>>().len() != indices.len()
                {
                    return Err(bad("Select bounded distinct donor annotation ordinals"));
                }
            }
        }
        bounded_json::preflight(self, MAX_REQUEST_BYTES)
    }
    /// Resolve actual retained ownership before copying a donor into a request.
    pub fn from_donor(
        workspace: &Workspace,
        event: &EventBinding,
        transform: Transform,
        indices: &[usize],
    ) -> Result<Self, Error> {
        workspace.validate()?;
        run::event_shape(event)?;
        let molecule = run::resolve(workspace, event)?;
        collagene_model::bounds::native_frame_preflight(
            std::iter::once(molecule),
            MAX_REQUEST_BYTES,
        )?;
        if indices.is_empty()
            || indices.len() > MAX_SELECTIONS
            || indices.iter().any(|i| *i >= molecule.features.len())
            || indices.iter().collect::<BTreeSet<_>>().len() != indices.len()
            || molecule.sequence.alphabet() != collagene_model::Alphabet::Dna
            || transform.origin >= molecule.sequence.len()
            || (molecule.topology == collagene_model::Topology::Linear && transform.origin != 0)
        {
            return Err(bad("Select bounded donor features"));
        }
        // Include request/container headroom before the model copy.
        bounded_json::preflight(molecule, MAX_REQUEST_BYTES - 64 * 1024)?;
        let source = Self::Donor {
            molecule: molecule.clone(),
            event: event.clone(),
            transform,
            indices: indices.to_vec(),
        };
        source.preflight()?;
        Ok(source)
    }
    /// Genuine event ownership admission; standalone source DTOs cannot grant it.
    pub fn validate_workspace(&self, workspace: &Workspace) -> Result<(), Error> {
        self.preflight()?;
        if let Self::Donor {
            molecule, event, ..
        } = self
        {
            workspace.validate()?;
            if run::resolve(workspace, event)? != molecule {
                return Err(bad("Donor differs from its exact retained workspace event"));
            }
        }
        Ok(())
    }
}
impl Request {
    pub fn preflight(&self) -> Result<(), Error> {
        if self.schema_version != 1 {
            return Err(Error::Unsupported(
                "Unknown annotation merge request schema".into(),
            ));
        }
        if !hex(&self.target_revision)
            || self.selections.is_empty()
            || self.selections.len() > MAX_SELECTIONS
        {
            return Err(bad("Merge requires exact target pin and 1..1000 decisions"));
        }
        if !(1..=library::MAX_SCAN_COMPARISONS).contains(&self.limits.max_comparisons)
            || !(1..=100_000).contains(&self.limits.max_relations)
            || !(4096..=MAX_REPORT_BYTES).contains(&self.limits.max_output_bytes)
        {
            return Err(Error::TooLarge);
        }
        self.source.preflight()?;
        let mut ids = BTreeSet::new();
        for selection in &self.selections {
            if !text(&selection.proposal_id, 128) || !ids.insert(&selection.proposal_id) {
                return Err(bad(
                    "Merge proposal identifiers must be distinct and bounded",
                ));
            }
            if let Decision::AddProvenance {
                target_index,
                annotation_revision,
            } = &selection.decision
            {
                if *target_index >= 100_000 || !hex(annotation_revision) {
                    return Err(bad("Invalid annotation attribution target pin"));
                }
            }
        }
        bounded_json::preflight(self, MAX_REQUEST_BYTES)
    }
    pub fn validate(&self, target: &Molecule) -> Result<(), Error> {
        self.preflight()?;
        collagene_model::bounds::native_frame_preflight(std::iter::once(target), MAX_REPORT_BYTES)?;
        target.validate()?;
        if target.revision() != self.target_revision {
            return Err(bad("Annotation merge target revision changed"));
        }
        match &self.source {
            Source::Library { library, report }
                if report.source_revision != self.target_revision
                    || report.library_revision != library.revision()
                    || report.engine_version != library::ENGINE_VERSION =>
            {
                return Err(bad("Library merge source/library/engine pin changed"));
            }
            Source::Dataset { report } if report.reference_revision != self.target_revision => {
                return Err(bad("Dataset merge target pin changed"));
            }
            _ => {}
        }
        for feature in &target.features {
            Ledger::from_feature(feature)?;
        }
        Ok(())
    }
}
impl Preview {
    pub(super) fn identity(&self) -> Result<String, Error> {
        bounded_json::hash(&(
            &self.engine_version,
            &self.context_revision,
            &self.target_revision,
            self.source_complete,
            &self.candidates,
            self.comparisons,
            self.complete,
            &self.warnings,
        ))
    }
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if self.engine_version != ENGINE_VERSION
            || self.context_revision != expected
            || !hex(expected)
            || !hex(&self.target_revision)
            || self.candidates.is_empty()
            || self.candidates.len() > MAX_SELECTIONS
            || self.comparisons > library::MAX_SCAN_COMPARISONS
            || self.warnings.len() > 256
        {
            return Err(bad("Invalid annotation merge review identity/bounds"));
        }
        if self.complete
            != (self.source_complete && self.candidates.iter().all(|c| c.relations_complete))
        {
            return Err(bad("Merge completeness disagrees with retained scope"));
        }
        let mut budget = collagene_model::bounds::ModelBudget::annotation_edit();
        let mut ids = BTreeSet::new();
        let mut relations = 0usize;
        for candidate in &self.candidates {
            budget.features(std::slice::from_ref(&candidate.feature))?;
            Ledger::from_feature(&candidate.feature)?;
            if !text(&candidate.proposal_id, 128) || !ids.insert(&candidate.proposal_id) {
                return Err(bad("Invalid duplicate merge candidate"));
            }
            candidate.contribution.validate()?;
            if candidate.contribution.feature_revision
                != provenance::semantic_revision(&candidate.feature)?
            {
                return Err(bad(
                    "Merge contribution differs from its proposed annotation",
                ));
            }
            relations = relations
                .checked_add(candidate.relations.len())
                .filter(|n| *n <= 100_000)
                .ok_or(Error::TooLarge)?;
            let mut target_indices = BTreeSet::new();
            for relation in &candidate.relations {
                if !target_indices.insert(relation.target_index) {
                    return Err(bad("Repeated target relationship ordinal"));
                }
                if relation.target_index >= 100_000 || !hex(&relation.annotation_revision) {
                    return Err(bad("Invalid target relationship pin"));
                }
            }
        }
        bounded_json::preflight(self, MAX_REPORT_BYTES)?;
        if self.report_revision != self.identity()? {
            return Err(bad("Annotation merge review contents changed"));
        }
        Ok(())
    }
    pub fn validate_request(&self, target: &Molecule, request: &Request) -> Result<(), Error> {
        self.validate_binding(&revision(target, request)?)?;
        if self.target_revision != request.target_revision
            || self.candidates.len() != request.selections.len()
            || self.comparisons > request.limits.max_comparisons
            || self
                .candidates
                .iter()
                .map(|c| c.relations.len())
                .sum::<usize>()
                > request.limits.max_relations
        {
            return Err(bad("Merge review differs from requested bounds"));
        }
        for (candidate, selection) in self.candidates.iter().zip(&request.selections) {
            if candidate.proposal_id != selection.proposal_id
                || candidate.decision != selection.decision
            {
                return Err(bad("Merge review decisions changed"));
            }
        }
        bounded_json::preflight(self, request.limits.max_output_bytes)
    }
    pub fn validate_for(&self, target: &Molecule, request: &Request) -> Result<(), Error> {
        self.validate_request(target, request)?;
        if self != &preview(target, request)? {
            return Err(bad("Merge review differs from fresh source assessment"));
        }
        Ok(())
    }
}
impl Applied {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.review.validate_binding(expected)?;
        collagene_model::bounds::native_frame_preflight(
            std::iter::once(&self.molecule),
            MAX_REPORT_BYTES,
        )?;
        self.molecule.validate()?;
        let retained = self.receipt.retained_original_features;
        if self.receipt.context_revision != expected
            || self.receipt.product_revision != self.molecule.revision()
            || self.receipt.actions.len() != self.review.candidates.len()
            || retained > self.molecule.features.len()
            || self.receipt.added_feature_indices.len() > MAX_SELECTIONS
            || self.receipt.comparisons < self.review.comparisons
            || self.receipt.comparisons > library::MAX_SCAN_COMPARISONS
        {
            return Err(bad("Invalid applied annotation merge receipt"));
        }
        bounded_json::preflight(self, MAX_REPORT_BYTES)?;
        for feature in &self.molecule.features {
            Ledger::from_feature(feature)?;
        }
        if self
            .receipt
            .added_feature_indices
            .iter()
            .copied()
            .ne(retained..self.molecule.features.len())
        {
            return Err(bad(
                "Added annotations must be the exact new ordinal suffix",
            ));
        }
        let mut next_append = retained;
        for (action, candidate) in self.receipt.actions.iter().zip(&self.review.candidates) {
            candidate
                .feature
                .location
                .validate(self.molecule.sequence.len())?;
            if action.proposal_id != candidate.proposal_id || action.decision != candidate.decision
            {
                return Err(bad("Applied merge action differs from review"));
            }
            match &action.decision {
                Decision::Skip => {
                    if action.destination_index.is_some() || action.contribution_added {
                        return Err(bad(
                            "Skipped annotation cannot have a destination/contribution",
                        ));
                    }
                }
                Decision::Append => {
                    if !candidate.relations_complete
                        || action.destination_index != Some(next_append)
                        || !action.contribution_added
                    {
                        return Err(bad(
                            "Appended annotation has inconsistent ordinal or review scope",
                        ));
                    }
                    let mut expected_feature = candidate.feature.clone();
                    expected_feature.needs_review = true;
                    Ledger::add_to(&mut expected_feature, &candidate.contribution)?;
                    if self.molecule.features.get(next_append) != Some(&expected_feature) {
                        return Err(bad(
                            "Appended annotation differs from the exact proposal and contribution",
                        ));
                    }
                    next_append += 1;
                }
                Decision::AddProvenance { target_index, .. } => {
                    if !candidate.relations_complete
                        || *target_index >= retained
                        || action.destination_index != Some(*target_index)
                        || !Ledger::from_feature(&self.molecule.features[*target_index])?
                            .contains_contribution(&candidate.contribution)
                    {
                        return Err(bad(
                            "Annotation attribution differs from its exact retained target",
                        ));
                    }
                }
            }
        }
        if next_append != self.molecule.features.len() {
            return Err(bad("Unrecorded appended annotation"));
        }
        bounded_json::preflight(self, MAX_REPORT_BYTES)
    }
    /// Exact fresh source/output binding without rerunning a matcher or importer.
    pub fn validate_request(&self, target: &Molecule, request: &Request) -> Result<(), Error> {
        self.review.validate_request(target, request)?;
        self.validate_binding(&self.review.context_revision)?;
        if self.receipt.comparisons > request.limits.max_comparisons
            || self.receipt.retained_original_features != target.features.len()
            || self.molecule.sequence != target.sequence
            || self.molecule.topology != target.topology
            || self.molecule.name != target.name
            || self.molecule.original_genbank != target.original_genbank
        {
            return Err(bad(
                "Annotation merge changed source calls/topology or retained source binding",
            ));
        }
        for (index, original) in target.features.iter().enumerate() {
            let output = &self.molecule.features[index];
            let mut expected = Ledger::from_feature(original)?;
            let mut changed = false;
            for (action, candidate) in self.receipt.actions.iter().zip(&self.review.candidates) {
                if matches!(&action.decision, Decision::AddProvenance{target_index,..} if *target_index==index)
                {
                    changed = true;
                    let new = !expected.contains_contribution(&candidate.contribution);
                    if new {
                        expected.contributions.push(candidate.contribution.clone());
                    }
                    if new != action.contribution_added {
                        return Err(bad("Contribution idempotency receipt changed"));
                    }
                }
            }
            if changed {
                expected.validate()?;
                if output.kind != original.kind
                    || output.location != original.location
                    || output.needs_review != original.needs_review
                    || output
                        .qualifiers
                        .iter()
                        .filter(|(k, _)| k != QUALIFIER)
                        .ne(original.qualifiers.iter().filter(|(k, _)| k != QUALIFIER))
                    || Ledger::from_feature(output)? != expected
                {
                    return Err(bad(
                        "Attribution altered existing annotation meaning or prior contributions",
                    ));
                }
            } else if output != original {
                return Err(bad("Merge altered an unselected retained annotation"));
            }
        }
        let dataset = match &request.source {
            Source::Dataset { report }
                if self.receipt.actions.iter().any(|a| a.contribution_added) =>
            {
                Some((
                    format!("annotation_import:{}:original_base64", report.source_blake3),
                    format!("annotation_import:{}:format", report.source_blake3),
                    report,
                ))
            }
            _ => None,
        };
        for (key, value) in &target.metadata {
            if self.molecule.metadata.get(key) != Some(value) {
                return Err(bad("Merge changed existing document metadata"));
            }
        }
        for (key, value) in &self.molecule.metadata {
            if target.metadata.contains_key(key) {
                continue;
            }
            let valid = dataset.as_ref().is_some_and(|(raw, format, report)| {
                (key == raw && value == &report.original_base64)
                    || (key == format && value == &format!("{:?}", report.format))
            });
            if !valid {
                return Err(bad("Merge added unrecorded document metadata"));
            }
        }
        bounded_json::preflight(self, request.limits.max_output_bytes)
    }
}
