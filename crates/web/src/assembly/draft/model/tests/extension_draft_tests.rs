use super::*;
use crate::cloning::preparation::cloning_methods_model::OverlapExtensionPcrInput;

#[test]
fn six_independent_unfinished_oligos_and_policy_roundtrip() {
    let mut draft = Draft::default();
    let mut raw = OverlapExtensionPcrInput::default();
    raw.left.forward.annealing = "ag tc unfinished".into();
    raw.left.reverse.tail = "GGTCTGAAAA".into();
    raw.right.forward.tail = "GCTACGTATT".into();
    raw.right.reverse.name = "Right outer stock primer".into();
    raw.outer.forward.tail = " separate final assembly tail ".into();
    raw.outer.reverse.annealing = "TGACGATC".into();
    raw.overlap = "1e".into();
    raw.max_candidates = String::new();
    draft.methods.rows[1].preparation = OverlapPreparation::OverlapExtensionPcr;
    draft.methods.rows[1].overlap_extension_pcr = Some(raw);
    draft.sync_rows();
    assert_eq!(draft.version, CURRENT_VERSION);
    let bytes = encode_raw(&draft).unwrap();
    let restored = decode(&bytes).unwrap();
    assert!(restored == draft);
    assert!(
        crate::cloning::preparation::cloning_preparations::prepared_rows(&restored.methods, 2)
            .is_err()
    );

    let mut old = restored.clone();
    old.version = CURRENT_VERSION - 1;
    assert!(old.validate().is_err());
    // Changing preparation preserves the inactive editable fields.
    old.methods.rows[1].preparation = OverlapPreparation::WholeLinear;
    assert!(old.validate().is_err());
    old.sync_rows();
    assert_eq!(old.version, CURRENT_VERSION);
    assert!(decode(&encode_raw(&old).unwrap()).unwrap() == old);
    let mut oversized = old;
    oversized.methods.rows[1]
        .overlap_extension_pcr
        .as_mut()
        .unwrap()
        .outer
        .reverse
        .tail = "A".repeat(1001);
    assert!(encode_raw(&oversized).is_err());
}

#[test]
fn physical_homology_settings_roundtrip_unfinished_inactive_values() {
    let mut draft = Draft::default();
    assert!(draft.methods.homology.is_none());
    let mut form = crate::cloning::preparation::cloning_methods_model::HomologyForm::default();
    form.enabled = false;
    form.exact_homology = false;
    form.cuts = "20:unfinished, 3:".into();
    form.outer_left = "partial".into();
    draft.methods.homology = Some(form);
    draft.sync_rows();
    assert_eq!(draft.version, CURRENT_VERSION);
    let restored = decode(&encode_raw(&draft).unwrap()).unwrap();
    assert!(restored == draft);
    let mut downgraded = restored.clone();
    downgraded.version = CURRENT_VERSION - 1;
    assert!(downgraded.validate().is_err());

    draft.methods.homology.as_mut().unwrap().cuts = "x".repeat(1025);
    assert!(encode_raw(&draft).is_err());
}
#[test]
fn three_prime_settings_preserve_independent_five_prime_settings() {
    let mut draft = Draft::default();
    assert!(draft.methods.annealing.is_none());
    draft.methods.homology =
        Some(crate::cloning::preparation::cloning_methods_model::HomologyForm::default());
    draft.sync_rows();
    assert_eq!(draft.version, CURRENT_VERSION);
    let current = decode(&encode_raw(&draft).unwrap()).unwrap();
    assert!(current.methods.annealing.is_none());
    assert!(current.methods.homology.as_ref().unwrap().enabled);
    let mut form = crate::cloning::preparation::cloning_methods_model::AnnealingForm::default();
    form.enabled = false;
    form.cuts = "0:unfinished, 20:".into();
    form.max_candidates = "partial".into();
    draft.methods.annealing = Some(form);
    draft.sync_rows();
    assert_eq!(draft.version, CURRENT_VERSION);
    let restored = decode(&encode_raw(&draft).unwrap()).unwrap();
    assert!(restored == draft);
    assert!(restored.methods.homology.as_ref().unwrap().enabled);
    assert!(!restored.methods.annealing.as_ref().unwrap().enabled);
    let mut downgraded = restored.clone();
    downgraded.version = CURRENT_VERSION - 1;
    assert!(downgraded.validate().is_err());

    draft.methods.annealing.as_mut().unwrap().cuts = "x".repeat(1025);
    assert!(encode_raw(&draft).is_err());
}
#[test]
fn extension_catalog_selector_preserves_unfinished_raw6() {
    let mut draft = Draft::default();
    draft.methods.rows[1].preparation = OverlapPreparation::OverlapExtensionPcrCatalogDigest;
    draft.sync_rows();
    assert_eq!(draft.version, CURRENT_VERSION);
    assert!(draft.methods.rows[1].overlap_extension_pcr.is_none());
    let mut raw = OverlapExtensionPcrInput::default();
    raw.outer.forward.tail = "CCGAATTC".into();
    raw.right.forward.annealing = "unfinished".into();
    draft.methods.rows[1].overlap_extension_pcr = Some(raw);
    draft.sync_rows();
    let restored = decode(&encode_raw(&draft).unwrap()).unwrap();
    assert!(restored == draft);

    let mut downgraded = draft.clone();
    downgraded.version = CURRENT_VERSION - 1;
    assert!(downgraded.validate().is_err());
    downgraded.methods.rows[1].preparation = OverlapPreparation::OverlapExtensionPcr;
    downgraded.version = CURRENT_VERSION;
    assert!(downgraded.validate().is_ok());
    assert!(decode(&encode_raw(&downgraded).unwrap()).unwrap() == downgraded);
}
#[test]
fn selected_empty_extension_preparation_roundtrips_explicit_null() {
    let mut draft = Draft::default();
    draft.methods.rows[0].preparation = OverlapPreparation::OverlapExtensionPcr;
    assert!(draft.validate().is_ok());
    draft.sync_rows();
    assert_eq!(draft.version, CURRENT_VERSION);
    assert!(draft.methods.rows[0].overlap_extension_pcr.is_none());
    assert!(decode(&encode_raw(&draft).unwrap()).unwrap() == draft);
}
