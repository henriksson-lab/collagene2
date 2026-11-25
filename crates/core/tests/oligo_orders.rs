use collagene_core::{
    Error,
    analysis::Primer,
    oligo_formats::{self, Delimiter, TableExportMode},
    oligo_orders,
};
#[test]
fn named_single_and_pair_orders_preserve_literal_rows_tails_and_spreadsheet_policy() {
    let one = Primer::new("=One, α", "ACGTR", "NN").unwrap();
    let single = oligo_orders::from_primers(&[one.clone()]).unwrap();
    assert_eq!(single.records.len(), 1);
    assert_eq!(single.records[0].name, "=One, α");
    assert_eq!(
        single.records[0].raw_ordered_sequence.as_deref(),
        Some("NNACGTR")
    );
    assert_eq!(
        single.records[0].raw_annealing_sequence.as_deref(),
        Some("ACGTR")
    );
    assert_eq!(single.records[0].raw_tail.as_deref(), Some("NN"));
    assert!(
        single.records[0]
            .diagnostics
            .iter()
            .any(|d| d.contains("IUPAC"))
    );
    let second = Primer::new("=One, α", "TCGA", "GC").unwrap();
    let pair = oligo_orders::from_primers(&[one, second]).unwrap();
    assert_eq!(pair.records.len(), 2);
    assert_ne!(pair.records[0].row_id, pair.records[1].row_id);
    assert_eq!(
        pair.records[1].raw_ordered_sequence.as_deref(),
        Some("GCTCGA")
    );
    pair.validate().unwrap();
    for delimiter in [Delimiter::Comma, Delimiter::Tab] {
        let safe = oligo_formats::export_oligos(&pair, delimiter, TableExportMode::SpreadsheetSafe)
            .unwrap();
        assert!(safe.text.contains("'=One, α"));
        assert!(safe.text.contains("NNACGTR"));
    }
    assert_eq!(pair.records[0].name, "=One, α");
    assert!(
        pair.records
            .iter()
            .all(|r| r.modification.is_none() && r.quantity.is_none())
    );
}
#[test]
fn long_order_evidence_and_existing_field_bounds_do_not_fabricate_reusable_primers() {
    let long = Primer::new("long order", "ACGA", &"A".repeat(250)).unwrap();
    let table = oligo_orders::from_primers(&[long]).unwrap();
    assert!(table.records[0].primer.is_none());
    assert_eq!(
        table.records[0]
            .raw_ordered_sequence
            .as_ref()
            .unwrap()
            .len(),
        254
    );
    table.validate().unwrap();
    assert!(oligo_orders::from_primers(&[]).is_err());
    assert!(oligo_orders::from_primers(&[Primer::new("\ninvalid", "ACGA", "").unwrap()]).is_err());
    assert!(matches!(
        oligo_orders::from_primers(&[Primer::new("oversize", "ACGA", &"T".repeat(8192)).unwrap()]),
        Err(Error::TooLarge)
    ));
}
