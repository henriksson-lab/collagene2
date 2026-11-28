use super::*;
#[test]
fn current_schema_requires_all_fields_including_explicit_nulls_and_defaults() {
    let draft = Draft::default();
    let value = serde_json::to_value(&draft).unwrap();
    assert!(decode(&serde_json::to_vec(&value).unwrap()).unwrap() == draft);
    for pointer in [
        "/junction_edits",
        "/substitutions",
        "/edits",
        "/alternatives",
        "/slots/0/selector",
        "/slots/0/captured",
        "/slots/0/digest_geometry",
        "/slots/0/anneal_bottom",
        "/methods/tailing_profile",
        "/methods/vector/action",
        "/methods/golden_gate",
        "/methods/golden_gate_pool",
        "/methods/multisite",
        "/methods/source_roles",
        "/methods/homology",
        "/methods/annealing",
        "/slots/0/material_origin",
        "/methods/rows/0/template_strand",
        "/methods/rows/0/run_origin",
        "/methods/rows/0/source_generation",
        "/methods/rows/0/overlap_extension_pcr",
        "/methods/rows/0/catalog_selection",
        "/methods/rows/0/mismatch_max_mismatches",
        "/methods/rows/0/anneal_options",
        "/methods/rows/0/anneal_candidate",
        "/methods/rows/0/anneal_bottom_document",
        "/methods/rows/0/anneal_source_context",
    ] {
        let mut missing = value.clone();
        let (parent, field) = pointer.rsplit_once('/').unwrap();
        missing
            .pointer_mut(parent)
            .unwrap()
            .as_object_mut()
            .unwrap()
            .remove(field);
        assert!(
            decode(&serde_json::to_vec(&missing).unwrap()).is_err(),
            "Missing field accepted: {pointer}"
        );
    }
    for version in [CURRENT_VERSION - 1, CURRENT_VERSION + 1] {
        let mut unsupported = value.clone();
        unsupported["version"] = version.into();
        assert!(decode(&serde_json::to_vec(&unsupported).unwrap()).is_err());
    }
    let mut unknown = value;
    unknown["methods"]["obsolete_control"] = true.into();
    assert!(decode(&serde_json::to_vec(&unknown).unwrap()).is_err());
}
