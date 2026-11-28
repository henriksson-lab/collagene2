use super::*;
#[test]
fn search_limits_and_unfinished_numeric_text_roundtrip_without_changing_intent() {
    let mut draft = Draft::default();
    let row = &mut draft.methods.rows[1];
    row.preparation = OverlapPreparation::MismatchPcrCatalogDigest;
    row.forward = "AGACCGAT".into();
    row.reverse_primer = "TCCAGGTC".into();
    row.mismatch_max_products = "1".into();
    row.mismatch_max_work = "1e".into();
    let restored = decode(&encode_raw(&draft).unwrap()).unwrap();
    assert!(restored == draft);
    assert!(restored.methods.rows[1].mismatch_options().is_err());
    draft.methods.rows[1].mismatch_max_work = "1000".into();
    let restored = decode(&encode_raw(&draft).unwrap()).unwrap();
    let options = restored.methods.rows[1].mismatch_options().unwrap();
    assert_eq!(options.max_products, 1);
    assert_eq!(options.binding.max_work, 1000);
    assert!(
        restored.slots == draft.slots
            && restored.junction_edits == draft.junction_edits
            && restored.substitutions == draft.substitutions
    );
    draft.methods.rows[1].preparation = OverlapPreparation::WholeLinear;
    draft.methods.rows[1].mismatch_max_pairs = "x".repeat(21);
    assert!(draft.validate().is_err());
}
