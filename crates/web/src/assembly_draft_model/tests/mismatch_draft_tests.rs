use super::*;
#[test]
fn unfinished_mismatch_primers_and_policy_are_persisted_without_claiming_valid_preparation() {
    let mut draft = Draft::default();
    draft.slots.truncate(1);
    draft.sync_rows();
    let row = &mut draft.methods.rows[0];
    row.preparation = OverlapPreparation::MismatchPcr;
    row.forward_name = "Unfinished mutation F".into();
    row.forward = String::new();
    row.reverse_primer = "GTTACCCC".into();
    row.mismatch_max_mismatches = "2e".into();
    row.mismatch_three_prime_exact = String::new();
    let restored = decode(&encode_raw(&draft).unwrap()).unwrap();
    assert!(restored == draft);
    assert!(crate::cloning_preparations::prepared_rows(&restored.methods, 1).is_err());
    draft.methods.rows[0].forward = "ACGACCCC".into();
    // Even with complete oligos, unfinished policy text remains an editable
    // draft while fresh preparation compilation reports the local error.
    assert!(decode(&encode_raw(&draft).unwrap()).unwrap() == draft);
    let error = crate::cloning_preparations::prepared_rows(&draft.methods, 1)
        .err()
        .unwrap();
    assert!(error.contains("allowed mismatches"));
}
#[test]
fn mismatch_policy_roundtrip_preserves_explicit_primers_and_unfinished_values() {
    let mut draft = Draft::default();
    let row = &mut draft.methods.rows[0];
    row.preparation = OverlapPreparation::MismatchPcr;
    row.forward = "ACGACCCC".into();
    row.reverse_primer = "GTTACCCC".into();
    row.mismatch_max_mismatches = "2".into();
    row.mismatch_three_prime_exact = "3".into();
    let restored = decode(&encode_raw(&draft).unwrap()).unwrap();
    assert!(restored == draft);
    draft.methods.rows[0].mismatch_three_prime_exact = "9".into();
    assert!(draft.validate().is_ok());
    assert!(draft.methods.rows[0].mismatch_options().is_err());
    assert!(decode(&encode_raw(&draft).unwrap()).unwrap() == draft);
    draft.methods.rows[0].preparation = OverlapPreparation::ExactPcr;
    draft.methods.rows[0].mismatch_max_mismatches = "x".repeat(21);
    assert!(draft.validate().is_err());
}
