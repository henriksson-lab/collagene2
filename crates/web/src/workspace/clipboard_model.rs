//! Clipboard text/edit descriptions scoped to the exact displayed molecule.
//! This module performs no document mutation, browser access or history writes.
use collagene_model::{Molecule, Sequence};

#[derive(Debug, PartialEq, Eq)]
pub struct Edit {
    pub start: usize,
    pub end: usize,
    pub replacement: String,
    pub caret: usize,
}
fn reference_range(
    molecule: &Molecule,
    span: Option<(usize, usize)>,
    parts: &[(usize, usize)],
    allow_caret: bool,
) -> Result<(usize, usize), String> {
    if parts.len() > 1 {
        return Err(
            "Choose one contiguous reference interval; multiple highlighted parts are ambiguous."
                .into(),
        );
    }
    let (start, end) = span.ok_or("Select a reference interval or place a caret first.")?;
    if parts.first().is_some_and(|part| *part != (start, end)) {
        return Err(
            "The highlighted part no longer matches the selected reference interval.".into(),
        );
    }
    if start > end || end > molecule.sequence.len() || (!allow_caret && start == end) {
        return Err("The selected reference interval is empty or no longer available.".into());
    }
    Ok((start, end))
}
fn editable(read_only: bool, locked: bool) -> Result<(), String> {
    if read_only {
        return Err("This displayed reference is read-only.".into());
    }
    if locked {
        return Err("Wait for the document update to finish.".into());
    }
    Ok(())
}
pub fn copy_text(
    molecule: &Molecule,
    feature: Option<usize>,
    span: Option<(usize, usize)>,
    parts: &[(usize, usize)],
) -> Result<String, String> {
    if let Some(index) = feature {
        return molecule
            .features
            .get(index)
            .ok_or("Selected annotation changed")?
            .location
            .extract(&molecule.sequence)
            .map_err(|error| error.to_string());
    }
    let (start, end) = reference_range(molecule, span, parts, false)?;
    molecule
        .sequence
        .as_str()
        .get(start..end)
        .map(str::to_owned)
        .ok_or_else(|| "Selected bases are no longer available.".into())
}
pub fn cut(
    molecule: &Molecule,
    feature: Option<usize>,
    span: Option<(usize, usize)>,
    parts: &[(usize, usize)],
    read_only: bool,
    locked: bool,
) -> Result<Edit, String> {
    editable(read_only, locked)?;
    if feature.is_some() {
        return Err("Select reference bases explicitly before cutting an annotation.".into());
    }
    let (start, end) = reference_range(molecule, span, parts, false)?;
    Ok(Edit {
        start,
        end,
        replacement: String::new(),
        caret: start,
    })
}
pub fn paste(
    molecule: &Molecule,
    feature: Option<usize>,
    span: Option<(usize, usize)>,
    parts: &[(usize, usize)],
    read_only: bool,
    locked: bool,
    text: &str,
) -> Result<Edit, String> {
    editable(read_only, locked)?;
    if feature.is_some() {
        return Err(
            "Select a reference interval or caret explicitly before replacing an annotation."
                .into(),
        );
    }
    let (start, end) = reference_range(molecule, span, parts, true)?;
    if text.len() > collagene_core::formats::sequence::MAX_IMPORT_BYTES {
        return Err("Clipboard text exceeds the import size limit.".into());
    }
    if text.trim().is_empty() {
        return Err("Clipboard text contains no sequence letters.".into());
    }
    let normalized =
        Sequence::new(molecule.sequence.alphabet(), text).map_err(|error| error.to_string())?;
    molecule
        .sequence
        .len()
        .checked_sub(end - start)
        .and_then(|retained| retained.checked_add(normalized.len()))
        .filter(|length| *length <= collagene_model::MAX_SEQUENCE_LEN)
        .ok_or("Pasted sequence would exceed the core sequence-length bound.")?;
    let caret = start
        .checked_add(normalized.len())
        .ok_or("Paste caret overflow.")?;
    Ok(Edit {
        start,
        end,
        replacement: normalized.as_str().into(),
        caret,
    })
}

pub fn insert(
    molecule: &Molecule,
    feature: Option<usize>,
    span: Option<(usize, usize)>,
    parts: &[(usize, usize)],
    read_only: bool,
    locked: bool,
    text: &str,
) -> Result<Edit, String> {
    editable(read_only, locked)?;
    if feature.is_some() {
        return Err("Choose an explicit reference interval or caret before inserting.".into());
    }
    let (start, _) = reference_range(molecule, span, parts, true)?;
    paste(
        molecule,
        None,
        Some((start, start)),
        &[],
        read_only,
        locked,
        text,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_model::{Alphabet, Feature, Location, Topology};
    fn dna(letters: &str) -> Molecule {
        Molecule::dna("Displayed reference", letters, Topology::Linear).unwrap()
    }
    #[test]
    fn copy_uses_the_displayed_source_or_product_and_never_a_global_document() {
        let source = dna("AACCGGTT");
        let product = dna("TTAGCAAA");
        let before = (source.clone(), product.clone());
        assert_eq!(
            copy_text(&source, None, Some((1, 5)), &[(1, 5)]).unwrap(),
            "ACCG"
        );
        assert_eq!(
            copy_text(&product, None, Some((1, 5)), &[(1, 5)]).unwrap(),
            "TAGC"
        );
        assert!(copy_text(&source, None, Some((0, 2)), &[(0, 2), (6, 8)]).is_err());
        assert!(copy_text(&source, None, Some((0, 2)), &[(2, 4)]).is_err());
        assert!(copy_text(&source, None, None, &[(0, 2)]).is_err());
        assert!(copy_text(&source, None, Some((3, 3)), &[]).is_err());
        assert_eq!((source, product), before);
    }
    #[test]
    fn annotation_copy_uses_biological_complement_join_order_instead_of_first_highlight() {
        let mut source = dna("AAGCTTAC");
        source.features.push(Feature {
            kind: "CDS".into(),
            location: Location::Complement(Box::new(Location::Join(vec![
                Location::range(0, 3),
                Location::range(5, 8),
            ]))),
            qualifiers: vec![],
            needs_review: false,
        });
        // AAG + TAC -> AAGTAC -> biological reverse complement GTA CTT.
        assert_eq!(
            copy_text(&source, Some(0), Some((0, 3)), &[(0, 3), (5, 8)]).unwrap(),
            "GTACTT"
        );
        assert!(copy_text(&source, Some(1), Some((0, 3)), &[]).is_err());
    }
    #[test]
    fn readonly_locked_features_compound_and_stale_ranges_cannot_be_edited() {
        let source = dna("AACCGGTT");
        let before = source.clone();
        assert!(
            paste(&source, None, None, &[], true, false, "Ω")
                .unwrap_err()
                .contains("read-only")
        );
        assert!(
            cut(&source, Some(999), None, &[], false, true)
                .unwrap_err()
                .contains("finish")
        );
        for (feature, span, parts) in [
            (Some(0), Some((1, 3)), vec![]),
            (None, Some((0, 2)), vec![(0, 2), (6, 8)]),
            (None, Some((0, 2)), vec![(3, 5)]),
            (None, Some((0, 9)), vec![]),
            (None, Some((5, 3)), vec![]),
            (None, None, vec![]),
        ] {
            assert!(cut(&source, feature, span, &parts, false, false).is_err());
            assert!(paste(&source, feature, span, &parts, false, false, "AC").is_err());
        }
        assert_eq!(source, before);
    }
    #[test]
    fn paste_normalizes_using_the_displayed_alphabet_and_cut_reports_the_exact_reference_edit() {
        let source = dna("AACCGGTT");
        assert_eq!(
            paste(
                &source,
                None,
                Some((2, 5)),
                &[(2, 5)],
                false,
                false,
                " a\nc g "
            )
            .unwrap(),
            Edit {
                start: 2,
                end: 5,
                replacement: "ACG".into(),
                caret: 5
            }
        );
        assert_eq!(
            paste(&source, None, Some((8, 8)), &[], false, false, "tt").unwrap(),
            Edit {
                start: 8,
                end: 8,
                replacement: "TT".into(),
                caret: 10
            }
        );
        assert_eq!(
            cut(&source, None, Some((1, 6)), &[], false, false).unwrap(),
            Edit {
                start: 1,
                end: 6,
                replacement: String::new(),
                caret: 1
            }
        );
        assert!(cut(&source, None, Some((2, 2)), &[], false, false).is_err());
        assert_eq!(
            insert(
                &source,
                None,
                Some((2, 6)),
                &[(2, 6)],
                false,
                false,
                " g t "
            )
            .unwrap(),
            Edit {
                start: 2,
                end: 2,
                replacement: "GT".into(),
                caret: 4
            }
        );
        assert!(
            insert(
                &source,
                None,
                Some((0, 2)),
                &[(0, 2), (6, 8)],
                false,
                false,
                "A"
            )
            .is_err()
        );
        for text in ["", " \t\n", "Ω", "AU", ">fasta\nAC"] {
            assert!(paste(&source, None, Some((1, 1)), &[], false, false, text).is_err());
        }
        let mut rna = dna("ACGT");
        rna.sequence = Sequence::new(Alphabet::Rna, "ACGU").unwrap();
        assert_eq!(
            paste(&rna, None, Some((0, 1)), &[], false, false, " ug a ")
                .unwrap()
                .replacement,
            "UGA"
        );
        assert!(paste(&rna, None, Some((0, 0)), &[], false, false, "T").is_err());
        let mut protein = dna("ACGT");
        protein.sequence = Sequence::new(Alphabet::Protein, "MEKW").unwrap();
        let edit = paste(&protein, None, Some((1, 3)), &[], false, false, " f\ny* ").unwrap();
        assert_eq!(edit.replacement, "FY*");
        assert_eq!(edit.caret, 4);
        let full = dna(&"A".repeat(collagene_model::MAX_SEQUENCE_LEN));
        assert!(paste(&full, None, Some((0, 0)), &[], false, false, "G").is_err());
        assert!(
            paste(
                &source,
                None,
                Some((0, 0)),
                &[],
                false,
                false,
                &" ".repeat(collagene_core::formats::sequence::MAX_IMPORT_BYTES + 1)
            )
            .unwrap_err()
            .contains("import size")
        );
    }
}
