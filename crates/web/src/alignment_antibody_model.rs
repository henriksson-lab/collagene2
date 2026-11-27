//! Immutable alignment/row/epoch guards; no synthetic document/event identity.
use collagene_core::{
    alignment_antibody as a, alignment_formats::AlignmentDocument,
    antibody_numbering::domain_candidates as d,
};
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Capture {
    pub epoch: u64,
    pub source: a::Source,
    pub options: d::Options,
    pub limits: d::Limits,
    pub selection: collagene_core::antibody_numbering::domain_tracks::Selection,
    pub revision: String,
}
impl Capture {
    pub fn prepare(
        epoch: u64,
        doc: &AlignmentDocument,
        source: a::Source,
        options: d::Options,
        limits: d::Limits,
        selection: collagene_core::antibody_numbering::domain_tracks::Selection,
    ) -> Result<Self, String> {
        let revision = a::request_revision(doc, &source, &options, &limits, &selection)
            .map_err(|e| e.to_string())?;
        Ok(Self {
            epoch,
            source,
            options,
            limits,
            selection,
            revision,
        })
    }
    pub fn guard(
        &self,
        epoch: u64,
        source: &a::Source,
        options: &d::Options,
        limits: &d::Limits,
        selection: &collagene_core::antibody_numbering::domain_tracks::Selection,
    ) -> Result<(), String> {
        if self.epoch != epoch
            || self.source != *source
            || self.options != *options
            || self.limits != *limits
            || self.selection != *selection
        {
            return Err(
                "Alignment, focused row, workspace or analysis settings changed; analyze again."
                    .into(),
            );
        }
        Ok(())
    }
    pub fn validate_report(&self, doc: &AlignmentDocument, r: &a::Report) -> Result<(), String> {
        r.validate_binding(&self.revision)
            .map_err(|e| e.to_string())?;
        r.validate_source(
            doc,
            &self.source,
            &self.options,
            &self.limits,
            &self.selection,
        )
        .map_err(|e| e.to_string())
    }
}
pub fn description(report: &a::Report, candidate: usize, column: usize) -> String {
    if column >= report.column_count {
        return "Alignment column is outside this report.".into();
    }
    let Ok(index) = report.residue_columns.binary_search(&column) else {
        return format!(
            "Alignment column {} is a gap in the analyzed row; no residue or numbering label.",
            column + 1
        );
    };
    let Some(row) = report.at_column(candidate, column) else {
        return format!(
            "Alignment column {} · row residue {} has no supported assignment in this candidate; not proof of an outside-domain residue.",
            column + 1,
            index + 1
        );
    };
    let selection = &report.analysis.selection;
    let region = row
        .region
        .as_ref()
        .map(crate::antibody_overlay::region)
        .unwrap_or("unresolved region");
    match &row.assignment {
        Some(a) => format!(
            "Alignment column {} · row residue {} {} · conditional {} {} · {} {} · {}",
            column + 1,
            index + 1,
            row.residue,
            crate::antibody_overlay::scheme_label(&selection.numbering_scheme),
            a.position,
            crate::antibody_overlay::definition_label(&selection.region_definition),
            region,
            crate::antibody_overlay::origin(&a.evidence)
        ),
        None => format!(
            "Alignment column {} · row residue {} {} · unresolved numbering · {} {}",
            column + 1,
            index + 1,
            row.residue,
            crate::antibody_overlay::definition_label(&selection.region_definition),
            region
        ),
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::{
        Alphabet,
        alignment_formats::{AlignmentFormat, import_alignment},
    };
    #[test]
    fn focused_columns_keep_row_scope_but_equal_rows_foreign_epochs_and_stale_settings_do_not() {
        let seq = "EVQLLESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCARDYGMDVWGQGTLVTVSS";
        let doc = import_alignment(
            format!(">a\n--{seq}--\n>b\n--{seq}--\n").as_bytes(),
            AlignmentFormat::AlignedFasta,
            Alphabet::Protein,
            "Rows",
        )
        .unwrap();
        let source = a::source(&doc, &doc.rows[0].id).unwrap();
        let other = a::source(&doc, &doc.rows[1].id).unwrap();
        let (o, l) = crate::antibody_model::preset(300, 1, 60, 32);
        let c = Capture::prepare(
            7,
            &doc,
            source.clone(),
            o.clone(),
            l.clone(),
            test_selection(),
        )
        .unwrap();
        let r = a::analyze(&doc, &source, &o, &l, &test_selection(), &mut || false).unwrap();
        c.validate_report(&doc, &r).unwrap();
        let mut selection = test_selection();
        selection.numbering_scheme =
            collagene_core::antibody_numbering::domain_tracks::NumberingScheme::Kabat1991;
        assert!(c.guard(7, &source, &o, &l, &selection).is_err());
        let mut definition = test_selection();
        definition.region_definition =
            collagene_core::antibody_numbering::domain_tracks::RegionDefinition::ChothiaCanonical;
        assert!(c.guard(7, &source, &o, &l, &definition).is_err());
        c.guard(7, &source, &o, &l, &test_selection()).unwrap();
        assert!(c.guard(7, &other, &o, &l, &test_selection()).is_err());
        assert!(c.guard(8, &source, &o, &l, &test_selection()).is_err());
        let mut options = o;
        options.max_fr4_mismatches = 0;
        assert!(
            c.guard(7, &source, &options, &l, &test_selection())
                .is_err()
        );
        assert!(description(&r, 0, 0).contains("gap"));
        assert!(description(&r, 0, 2).contains("IMGT 1"));
        assert!(description(&r, 0, doc.column_count()).contains("outside this report"));
    }
}

#[cfg(test)]
fn test_selection() -> collagene_core::antibody_numbering::domain_tracks::Selection {
    collagene_core::antibody_numbering::domain_tracks::Selection {
        numbering_scheme: collagene_core::antibody_numbering::domain_tracks::NumberingScheme::Imgt,
        region_definition:
            collagene_core::antibody_numbering::domain_tracks::RegionDefinition::Imgt,
    }
}
