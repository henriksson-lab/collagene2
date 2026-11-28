//! Resource, immutable-source and documentary binding owners; no scientific replay.
use super::*;
use crate::{
    alignment::review,
    documents::workspace::Workspace,
    encoding::bounded_json,
    material::{lineage, run},
    sequence::translation,
};
use collagene_model::{Alphabet, Error, Molecule};
use std::collections::BTreeSet;
pub(super) fn bad(message: &str) -> Error {
    Error::Format(message.into())
}
fn hex(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || matches!(b, b'a'..=b'f'))
}
fn text(value: &str, max: usize) -> bool {
    !value.trim().is_empty() && value.len() <= max && !value.chars().any(char::is_control)
}
impl Request {
    pub(super) fn shape(&self) -> Result<(), Error> {
        self.limits.validate()?;
        if self.schema_version != 1 {
            return Err(Error::Unsupported(
                "Unsupported design/consensus request schema".into(),
            ));
        }
        run::event_shape(&self.design)?;
        run::event_shape(&self.consensus)?;
        if self.design == self.consensus {
            return Err(bad(
                "Design and supplied consensus require distinct exact retained events",
            ));
        }
        review::Options::Whole(self.comparison).validate()?;
        if self.cds.len() > MAX_CDS {
            return Err(Error::TooLarge);
        }
        if let ConsensusAxis::DeclaredSingleCircularTurn { rationale } = &self.consensus_axis {
            if !text(rationale, 2048) {
                return Err(bad("Declare a bounded circular-coordinate rationale"));
            }
        }
        let mut ids = BTreeSet::new();
        for selection in &self.cds {
            let assumption = &selection.design_cds;
            if !text(&selection.id, 96)
                || !ids.insert(&selection.id)
                || !hex(&assumption.source_revision)
                || !hex(&assumption.annotation_revision)
                || !text(&assumption.rationale, 2048)
            {
                return Err(bad(
                    "Invalid duplicate/unbounded CDS selection or annotation pin",
                ));
            }
            translation::CodeDescriptor::for_id(selection.genetic_code)?;
        }
        bounded_json::preflight(self, MAX_REQUEST_BYTES)
    }
    /// Exact event resolution for presentation; success is not scientific authority.
    /// Full validate admits the workspace and native source models separately.
    pub fn sources<'a>(
        &self,
        workspace: &'a Workspace,
    ) -> Result<(&'a Molecule, &'a Molecule), Error> {
        run::event_shape(&self.design)?;
        run::event_shape(&self.consensus)?;
        Ok((
            run::resolve(workspace, &self.design)?,
            run::resolve(workspace, &self.consensus)?,
        ))
    }
    pub fn validate(&self, workspace: &Workspace) -> Result<(), Error> {
        self.shape()?;
        workspace.validate()?;
        let (design, consensus) = self.sources(workspace)?;
        collagene_model::bounds::model_preflight(&[design, consensus])?;
        design.validate()?;
        consensus.validate()?;
        if design.sequence.alphabet() != Alphabet::Dna
            || consensus.sequence.alphabet() != Alphabet::Dna
        {
            return Err(Error::Unsupported(
                "Supplied whole-plasmid consensus assessment currently requires DNA sources".into(),
            ));
        }
        if design.sequence.len() > crate::alignment::MAX_ALIGNMENT_COLUMNS
            || consensus.sequence.len() > crate::alignment::MAX_ALIGNMENT_COLUMNS
        {
            return Err(Error::TooLarge);
        }
        self.comparison
            .alignment
            .validate_for_lengths(design.sequence.len(), consensus.sequence.len())
    }
}
pub(super) fn revision(workspace: &Workspace, request: &Request) -> Result<String, Error> {
    request.validate(workspace)?;
    bounded_json::hash(&(
        ENGINE_VERSION,
        lineage::workspace_fingerprint(workspace)?,
        request,
    ))
}
impl Evidence {
    fn identity(&self) -> Result<String, Error> {
        bounded_json::hash(&(
            &self.engine_version,
            &self.context_revision,
            &self.workspace_fingerprint,
            &self.request_revision,
            &self.pairwise_context_revision,
            &self.pairwise_report_revision,
            &self.request,
            &self.provenance,
            &self.cds,
            &self.work,
            self.coding_complete,
            self.selected_path_conditional,
            &self.statement,
        ))
    }
    pub(super) fn seal(&mut self, pairwise: &review::Report) -> Result<(), Error> {
        self.pairwise_context_revision = pairwise.context_revision.clone();
        self.pairwise_report_revision = bounded_json::hash(pairwise)?;
        self.evidence_revision = self.identity()?;
        Ok(())
    }
    /// Rebind retained semantic evidence to historical immutable events. This
    /// deliberately does not compare a pre-capture workspace fingerprint.
    pub fn validate_sources(
        &self,
        workspace: &Workspace,
        pairwise: &review::Report,
    ) -> Result<(), Error> {
        self.validate_for_pairwise(pairwise)?;
        self.request.validate(workspace)?;
        let (design, consensus) = self.request.sources(workspace)?;
        pairwise.validate_sources(design, consensus)?;
        let actual = [
            super::provenance::capture(SourceRole::Design, &self.request.design, design),
            super::provenance::capture(SourceRole::Consensus, &self.request.consensus, consensus),
        ];
        if self.provenance.as_slice() != actual.as_slice() {
            return Err(bad(
                "Retained source provenance differs from its immutable events",
            ));
        }
        for row in &self.cds {
            if row.design_translation.is_some() {
                let feature = design
                    .features
                    .get(row.design_cds.feature_index)
                    .ok_or_else(|| bad("Retained translated CDS feature is unavailable"))?;
                if row.design_cds.source_revision != self.request.design.revision
                    || crate::construct_qc::annotation_revision(feature)?
                        != row.design_cds.annotation_revision
                {
                    return Err(bad("Retained translated CDS annotation/source pin changed"));
                }
            }
        }
        Ok(())
    }
    /// Current bounded documentary/path binding, never a fresh semantic assessment.
    pub fn validate_for_pairwise(&self, pairwise: &review::Report) -> Result<(), Error> {
        self.request.shape()?;
        if self.engine_version != ENGINE_VERSION
            || !hex(&self.context_revision)
            || !hex(&self.workspace_fingerprint)
            || self.request_revision != bounded_json::hash(&self.request)?
            || self.context_revision
                != bounded_json::hash(&(
                    ENGINE_VERSION,
                    &self.workspace_fingerprint,
                    &self.request,
                ))?
            || self.statement != STATEMENT
            || !self.selected_path_conditional
            || self.cds.len() != self.request.cds.len()
            || self.provenance.len() != 2
            || self.work.mapping_work > self.request.limits.max_mapping_work
            || self.work.coding_calls > self.request.limits.max_coding_calls
            || self.work.segments > self.request.limits.max_segments
            || self.work.consequences > self.request.limits.max_consequences
        {
            return Err(bad(
                "Design/consensus evidence identity, scope or resource claims changed",
            ));
        }
        pairwise.validate_binding(&self.pairwise_context_revision)?;
        if pairwise.reference_revision != self.request.design.revision
            || pairwise.query_revision != self.request.consensus.revision
            || pairwise.options != review::Options::Whole(self.request.comparison)
        {
            return Err(bad(
                "Actual pairwise path differs from requested sources/options",
            ));
        }
        let mut segments = 0usize;
        let mut consequences = 0usize;
        for (row, selection) in self.cds.iter().zip(&self.request.cds) {
            if row.id != selection.id
                || row.design_cds != selection.design_cds
                || row.code != crate::construct_qc::CodePin::for_id(selection.genetic_code)?
                || row.register_segments.len() > 512
                || row.consequences.len() > 4096
                || row.findings.len() > 16
                || row
                    .findings
                    .iter()
                    .any(|finding| !text(&finding.message, 2048))
            {
                return Err(bad(
                    "Conditional CDS selection or retained evidence bounds changed",
                ));
            }
            segments = segments
                .checked_add(row.register_segments.len())
                .ok_or(Error::TooLarge)?;
            consequences = consequences
                .checked_add(row.consequences.len())
                .ok_or(Error::TooLarge)?;
            for segment in &row.register_segments {
                if segment.design_start >= segment.design_end
                    || segment.projected_start >= segment.projected_end
                    || segment.offset_mod3 > 2
                {
                    return Err(bad("Invalid conditional coding register segment"));
                }
            }
            for consequence in &row.consequences {
                if consequence.projected_codon.len() != 3
                    || consequence
                        .design_codon
                        .as_ref()
                        .is_some_and(|codon| codon.len() != 3)
                    || consequence.design_codon_index.is_some()
                        != consequence.design_codon.is_some()
                    || consequence.design_codon_index.is_some()
                        != consequence.design_amino_acid.is_some()
                    || consequence.design_coordinates.len() > 3
                    || consequence.consensus_coordinates.len() != 3
                    || !text(&consequence.message, 2048)
                {
                    return Err(bad("Invalid conditional codon evidence"));
                }
                for (span, length) in consequence
                    .design_coordinates
                    .iter()
                    .map(|span| (span, pairwise.reference_length))
                    .chain(
                        consequence
                            .consensus_coordinates
                            .iter()
                            .map(|span| (span, pairwise.query_length)),
                    )
                {
                    if span.start >= span.end || span.end > length {
                        return Err(bad("Codon coordinates outside exact retained sources"));
                    }
                }
            }
            for translation in row
                .design_translation
                .iter()
                .chain(&row.projected_translation)
            {
                if translation.protein.len() > self.request.limits.max_coding_calls / 3 + 1
                    || translation.options.genetic_code != row.code.descriptor.id
                    || translation.warnings.len() > 128
                    || translation
                        .warnings
                        .iter()
                        .any(|warning| warning.len() > 2048)
                {
                    return Err(Error::TooLarge);
                }
            }
        }
        if segments != self.work.segments
            || consequences != self.work.consequences
            || self.coding_complete
                != self.cds.iter().all(|row| {
                    row.mapping_complete && row.register_complete && row.consequences_complete
                })
        {
            return Err(bad("Coding completeness/work counters changed"));
        }
        for (index, provenance) in self.provenance.iter().enumerate() {
            let (event, role) = if index == 0 {
                (&self.request.design, SourceRole::Design)
            } else {
                (&self.request.consensus, SourceRole::Consensus)
            };
            if &provenance.event != event
                || provenance.role != role
                || provenance.topology
                    != if index == 0 {
                        pairwise.comparison.reference_topology
                    } else {
                        pairwise.comparison.query_topology
                    }
                || provenance.fields.len() > 8
                || provenance
                    .fields
                    .iter()
                    .any(|field| !text(&field.key, 64) || !text(&field.value, 2048))
                || provenance.original.byte_count > 4 * 1024 * 1024
                || !text(&provenance.original.message, 2048)
                || provenance
                    .original
                    .format
                    .as_ref()
                    .is_some_and(|format| !text(format, 32))
                || provenance
                    .original
                    .blake3
                    .as_ref()
                    .is_some_and(|hash| !hex(hash))
            {
                return Err(bad(
                    "Supplied-source provenance bounds or role binding changed",
                ));
            }
        }
        bounded_json::preflight(self, self.request.limits.max_report_bytes)?;
        if self.pairwise_report_revision != bounded_json::hash(pairwise)?
            || self.evidence_revision != self.identity()?
        {
            return Err(bad(
                "Conditional evidence no longer binds the actual retained pairwise content",
            ));
        }
        Ok(())
    }
}
impl Report {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if self.evidence.context_revision != expected {
            return Err(bad("Design/consensus context changed"));
        }
        self.evidence.validate_for_pairwise(&self.pairwise)?;
        bounded_json::preflight(self, self.evidence.request.limits.max_report_bytes)
    }
    /// Bounded exact-request/source binding; no alignment or coding replay.
    pub fn validate_request(&self, workspace: &Workspace, request: &Request) -> Result<(), Error> {
        let expected = revision(workspace, request)?;
        self.validate_binding(&expected)?;
        if &self.evidence.request != request
            || self.evidence.workspace_fingerprint != lineage::workspace_fingerprint(workspace)?
        {
            return Err(bad("Semantic request/workspace binding changed"));
        }
        let (design, consensus) = request.sources(workspace)?;
        self.pairwise.validate_sources(design, consensus)
    }
    /// Explicit heavy reassessment of imported evidence.
    pub fn validate_for(&self, workspace: &Workspace, request: &Request) -> Result<(), Error> {
        self.validate_request(workspace, request)?;
        if self.evidence.evidence_revision
            != super::assess(workspace, request)?
                .evidence
                .evidence_revision
        {
            return Err(bad(
                "Semantic evidence differs from fresh native assessment",
            ));
        }
        Ok(())
    }
}
