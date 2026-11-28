//! Independently expected lower-row references and coordinates; not Java output
//! execution or chemistry reconstructed from a serialized single strand.
use collagene_core::{
    chemistry::Duplex,
    collagene_xml_compat::{self as compat, ReferencePolicy as Policy},
    legacy_formats as legacy, protocol,
};
use collagene_model::{Location, Topology};
fn record(upper: &str, lower: &str, circular: bool, children: &str) -> String {
    format!(
        "<oneseq><sequence name=\"literal lower\" su=\"{upper}\" sl=\"{lower}\" circular=\"{circular}\" notes=\"literal &amp; retained\"/>{children}</oneseq>"
    )
}
fn annotation(from: i64, to: i64, direction: &str) -> String {
    format!(
        "<annotation name=\"footprint\" from=\"{from}\" to=\"{to}\" orientation=\"{direction}\" colr=\"1\" colg=\"2\" colb=\"3\"/>"
    )
}
#[test]
fn ordinary_profile_is_identical_and_mixed_lower_project_requires_explicit_review() {
    let ordinary = record("ACGTA", "TGCAT", false, "");
    assert!(!compat::needs_lower_only_review(ordinary.as_bytes()).unwrap());
    let old = legacy::import_collagene_xml(ordinary.as_bytes()).unwrap();
    for policy in [Policy::UpperOnly, Policy::AllowLowerOnly] {
        assert_eq!(
            serde_json::to_value(compat::import_with_policy(ordinary.as_bytes(), policy).unwrap())
                .unwrap(),
            serde_json::to_value(&old).unwrap()
        );
    }
    let lower = record("      ", "ACGTTA", false, "");
    let bytes = format!("<madgene>{ordinary}{lower}</madgene>").into_bytes();
    assert!(compat::needs_lower_only_review(&bytes).unwrap());
    assert!(compat::import_with_policy(&bytes, Policy::UpperOnly).is_err());
    let report = compat::import_with_policy(&bytes, Policy::AllowLowerOnly).unwrap();
    assert_eq!(report.molecules.len(), 2);
    assert_eq!(report.molecules[0].sequence.as_str(), "ACGTA");
    assert_eq!(report.molecules[1].sequence.as_str(), "ATTGCA");
    assert_eq!(
        report.molecules[1].metadata["description"],
        "literal & retained"
    );
    for (index, molecule) in report.molecules.iter().enumerate() {
        assert_eq!(molecule.metadata["legacy:record_index"], index.to_string());
        assert_eq!(legacy::original_legacy_bytes(molecule).unwrap(), bytes);
    }
    let source = &report.molecules[1];
    assert_eq!(source.metadata["legacy:upper_aligned"], "      ");
    assert_eq!(source.metadata["legacy:lower_aligned"], "ACGTTA");
    assert_eq!(source.metadata["legacy:reference_strand"], "lower-5prime");
    assert_eq!(source.metadata["genbank.molecule_type"], "ss-DNA");
    assert!(Duplex::from_molecule(source).is_err());
    let native = protocol::load_document(&protocol::save_document(source).unwrap()).unwrap();
    assert_eq!(&native, source);
    assert_eq!(legacy::original_legacy_bytes(&native).unwrap(), bytes);
}
#[test]
fn reversed_lower_reference_maps_only_literal_residue_footprints_and_keeps_oligos() {
    let children = format!(
        "{}{}{}<primer name=\"literal mismatching oligo\" sequence=\"AC\" target=\"1\" orientation=\"rev\"/>",
        annotation(0, 2, "fwd"),
        annotation(2, 5, "rev"),
        annotation(1, 3, "none")
    );
    let bytes = record("      ", "ACGTTA", false, &children);
    let report = compat::import_with_policy(bytes.as_bytes(), Policy::AllowLowerOnly).unwrap();
    let molecule = &report.molecules[0];
    assert_eq!(molecule.sequence.as_str(), "ATTGCA");
    assert_eq!(molecule.features[0].location.spans(), vec![(4, 6, true)]);
    assert_eq!(molecule.features[1].location.spans(), vec![(1, 4, false)]);
    assert_eq!(
        molecule.features[1]
            .location
            .extract(&molecule.sequence)
            .unwrap(),
        "TTG"
    );
    assert!(matches!(molecule.features[2].location, Location::Order(_)));
    assert!(
        molecule.features[2]
            .location
            .extract(&molecule.sequence)
            .is_err()
    );
    let primer = &molecule.features[3];
    assert_eq!(primer.location.spans(), vec![(3, 5, false)]);
    assert!(
        primer
            .qualifiers
            .iter()
            .any(|(key, value)| key == "oligo_sequence" && value.as_deref() == Some("AC"))
    );
    assert!(molecule.features.iter().all(|feature| feature.needs_review));
    assert_eq!(molecule.features.len(), 4);
    assert!(report.warnings.iter().any(|w| w.contains("No binding")));
}
#[test]
fn circular_wrap_is_reversed_in_order_and_gaps_or_ambiguous_ranges_stay_unresolved() {
    let bytes = record("      ", "ACGTTA", true, &annotation(4, 2, "rev"));
    let report = compat::import_with_policy(bytes.as_bytes(), Policy::AllowLowerOnly).unwrap();
    let molecule = &report.molecules[0];
    assert_eq!(molecule.topology, Topology::Circular);
    assert_eq!(
        molecule.features[0].location.spans(),
        vec![(4, 6, false), (0, 2, false)]
    );
    assert_eq!(
        molecule.features[0]
            .location
            .extract(&molecule.sequence)
            .unwrap(),
        "CAAT"
    );
    let children = format!(
        "{}{}{}{}",
        annotation(1, 3, "rev"),
        annotation(0, 2, "rev"),
        annotation(-1, 3, "rev"),
        annotation(2, 2, "fwd")
    );
    let bytes = record("      ", " AC G ", false, &children);
    let report = compat::import_with_policy(bytes.as_bytes(), Policy::AllowLowerOnly).unwrap();
    let molecule = &report.molecules[0];
    assert_eq!(molecule.sequence.as_str(), "GCA");
    assert_eq!(molecule.features[0].location.spans(), vec![(1, 3, false)]);
    assert_eq!(
        molecule.features[0]
            .location
            .extract(&molecule.sequence)
            .unwrap(),
        "CA"
    );
    assert!(
        molecule.features[1..]
            .iter()
            .all(|feature| matches!(feature.location, Location::Unresolved(_)))
    );
    assert!(Duplex::from_molecule(molecule).is_err());
}
#[test]
fn original_entities_long_attribute_spacing_and_unknown_workflow_are_retained_without_authority() {
    let bytes = format!(
        "<f><sequence name='é literal' su{}= {}'&#32;&#32;&#32;&#32;' sl='A&#67;G&#84;' circular='false' notes='literal &amp; retained'/><workflow method='unknown' chemistry='uninterpreted'/></f>",
        " ".repeat(300),
        " ".repeat(300)
    );
    let report = compat::import_with_policy(bytes.as_bytes(), Policy::AllowLowerOnly).unwrap();
    let molecule = &report.molecules[0];
    assert_eq!(molecule.name, "é literal");
    assert_eq!(molecule.sequence.as_str(), "TGCA");
    assert_eq!(
        legacy::original_legacy_bytes(molecule).unwrap(),
        bytes.as_bytes()
    );
    assert_eq!(
        molecule.metadata["source:legacy:chemistry_unmapped"],
        "true"
    );
    assert!(report.warnings.iter().any(|w| w.contains("workflow")));
    assert!(Duplex::from_molecule(molecule).is_err());
}
#[test]
fn malformed_or_over_budget_projects_fail_atomically_without_partial_record_recovery() {
    let lower = record("    ", "ACGT", false, "");
    for bytes in [
        format!("<madgene>{lower}<unknown/></madgene>"),
        format!("<madgene>{lower}<oneseq><sequence name='bad' su='ACGT' sl='TGCA' circular='maybe'/></oneseq></madgene>"),
        "<!DOCTYPE f [<!ENTITY row 'ACGT'>]><f><sequence name='bad' su='    ' sl='&row;' circular='false'/></f>".into(),
        record("    ","ACG!",false,""),
        format!("<madgene>{lower}<oneseq><sequence name='empty' su='    ' sl='    ' circular='false'/></oneseq></madgene>"),
        format!("<madgene>{lower}<oneseq><sequence name='bad' su='AC' sl='TG' circular='false'/><primer name='invalid' sequence='A!' target='1' orientation='fwd'/></oneseq></madgene>"),
    ]{assert!(compat::import_with_policy(bytes.as_bytes(),Policy::AllowLowerOnly).is_err());}
    let oversized = format!("<madgene>{}</madgene>", lower.repeat(1000));
    assert!(compat::import_with_policy(oversized.as_bytes(), Policy::AllowLowerOnly).is_err());
    let deep = format!(
        "<f><sequence name='lower' su='    ' sl='ACGT' circular='false'/>{}{}{}</f>",
        "<x>".repeat(65),
        "literal",
        "</x>".repeat(65)
    );
    assert!(compat::needs_lower_only_review(deep.as_bytes()).is_err());
}
