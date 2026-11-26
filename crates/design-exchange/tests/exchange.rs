use collagene_design_exchange::*;
use serde_json::{Value, json};
const HIERARCHY: &[u8] = include_bytes!("fixtures/hierarchy.jsonld");
const ALPHABETS: &[u8] = include_bytes!("fixtures/alphabets.jsonld");
fn selected() -> Selection {
    Selection {
        component: "https://example.org/lab/plasmid".into(),
        sequence: "https://example.org/lab/plasmid_sequence".into(),
        topology: None,
    }
}
fn bytes(v: &Value) -> Vec<u8> {
    serde_json::to_vec(v).unwrap()
}
fn object_mut<'a>(v: &'a mut Value, id: &str) -> &'a mut Value {
    v["@graph"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .find(|o| o["@id"] == id)
        .unwrap()
}
#[test]
fn real_component_sequence_hierarchy_and_constraints_are_retained_not_flattened() {
    let doc = import(HIERARCHY).unwrap();
    assert_eq!(doc.original_bytes(), HIERARCHY);
    let Kind::Component(c) = &doc.objects()[&selected().component].kind else {
        panic!()
    };
    assert_eq!(c.features.len(), 2);
    assert_eq!(c.constraints.len(), 1);
    let Kind::Feature(f) = &doc.objects()["https://example.org/lab/insert"].kind else {
        panic!()
    };
    assert!(
        matches!(&f.kind,FeatureKind::SubComponent{instance_of,source_locations,..}if instance_of=="https://example.org/lab/part"&&source_locations.len()==1)
    );
    assert!(doc.project(&selected(), ProjectionPolicy::Strict).is_err());
    let p = doc
        .project(&selected(), ProjectionPolicy::AllowOmissions)
        .unwrap();
    assert_eq!(p.molecule.sequence.as_str(), "AAAACCCCGGGGTTTT");
    assert_eq!(p.molecule.topology, Topology::Circular);
    assert_eq!(
        p.molecule.features[0]
            .location
            .extract(&p.molecule.sequence)
            .unwrap(),
        "CCCC"
    );
    assert_eq!(
        p.molecule.features[1]
            .location
            .extract(&p.molecule.sequence)
            .unwrap(),
        "TTTTAAAA"
    );
    assert!(!p.issues.is_empty());
    assert!(p.molecule.features.iter().all(|f| f.needs_review));
    assert!(!p.molecule.metadata.values().any(|v| v == "verified"));
    assert!(!p.molecule.metadata.contains_key("duplex_model"));
    let expanded = doc.expanded_bytes().unwrap();
    let again = import(&expanded).unwrap();
    assert_eq!(again.objects(), doc.objects());
    assert!(
        String::from_utf8(expanded)
            .unwrap()
            .contains("https://example.org/foreign/status")
    );
}
#[test]
fn explicit_rna_protein_and_case_normalization_never_guess_nucleic_alphabet() {
    let d = import(ALPHABETS).unwrap();
    let r = d
        .project(
            &Selection {
                component: "urn:example:rna".into(),
                sequence: "urn:example:rna_sequence".into(),
                topology: None,
            },
            ProjectionPolicy::Strict,
        )
        .unwrap();
    assert_eq!(r.molecule.sequence.alphabet(), Alphabet::Rna);
    assert_eq!(r.molecule.sequence.as_str(), "ACGU");
    assert_eq!(r.warnings.len(), 1);
    let p = d
        .project(
            &Selection {
                component: "urn:example:protein".into(),
                sequence: "urn:example:protein_sequence".into(),
                topology: None,
            },
            ProjectionPolicy::Strict,
        )
        .unwrap();
    assert_eq!(p.molecule.sequence.alphabet(), Alphabet::Protein);
    assert_eq!(p.molecule.sequence.as_str(), "MKWVTFISLL");
    let mut wrong: Value = serde_json::from_slice(ALPHABETS).unwrap();
    wrong[0][format!("{SBOL}type")] = json!([{"@id":DNA}]);
    assert!(
        import(&bytes(&wrong))
            .unwrap()
            .project(
                &Selection {
                    component: "urn:example:rna".into(),
                    sequence: "urn:example:rna_sequence".into(),
                    topology: Some(Topology::Linear)
                },
                ProjectionPolicy::AllowOmissions
            )
            .is_err()
    );
}
#[test]
fn unsupported_contexts_duplicate_properties_relative_ids_and_keyword_remapping_are_refused() {
    let mut v: Value = serde_json::from_slice(HIERARCHY).unwrap();
    v["@context"] = json!("https://attacker.invalid/context");
    assert!(import(&bytes(&v)).is_err());
    v["@context"] = json!({"sbol":"https://attacker.invalid/"});
    assert!(import(&bytes(&v)).is_err());
    v["@context"] = json!({"sbol":SBOL,"@vocab":SBOL});
    assert!(import(&bytes(&v)).is_err());
    assert!(import(br#"{"@graph":[],"@graph":[]}"#).is_err());
    let mut v: Value = serde_json::from_slice(HIERARCHY).unwrap();
    object_mut(&mut v, "https://example.org/lab/insert")["@id"] = json!("relative_part");
    assert!(import(&bytes(&v)).is_err());
    let mut v: Value = serde_json::from_slice(HIERARCHY).unwrap();
    object_mut(&mut v, "https://example.org/lab/plasmid")[format!("{SBOL}name")] =
        json!("aliased duplicate");
    assert!(import(&bytes(&v)).is_err());
    let mut v: Value = serde_json::from_slice(HIERARCHY).unwrap();
    object_mut(&mut v, "https://example.org/lab/insert_range")["sbol:start"] = json!({"@list":[1]});
    assert!(import(&bytes(&v)).is_err());
}
#[test]
fn wrong_sequence_ownership_cross_component_constraints_and_cyclic_hierarchies_do_not_pass() {
    let mut v: Value = serde_json::from_slice(HIERARCHY).unwrap();
    object_mut(&mut v, "https://example.org/lab/insert_range")["sbol:hasSequence"] =
        json!({"@id":"https://example.org/lab/part_sequence"});
    assert!(import(&bytes(&v)).is_err());
    let mut v: Value = serde_json::from_slice(HIERARCHY).unwrap();
    object_mut(&mut v, "https://example.org/lab/ordering")["sbol:object"] =
        json!({"@id":"https://example.org/lab/part"});
    assert!(import(&bytes(&v)).is_err());
    let mut v: Value = serde_json::from_slice(HIERARCHY).unwrap();
    object_mut(&mut v, "https://example.org/lab/part")["sbol:hasFeature"] =
        json!({"@id":"https://example.org/lab/recurse"});
    v["@graph"].as_array_mut().unwrap().push(json!({"@id":"https://example.org/lab/recurse","@type":"sbol:SubComponent","sbol:instanceOf":{"@id":"https://example.org/lab/plasmid"}}));
    assert!(import(&bytes(&v)).is_err());
    let mut v: Value = serde_json::from_slice(HIERARCHY).unwrap();
    object_mut(&mut v, "https://example.org/lab/insert_range")["sbol:end"] = json!(4);
    assert!(import(&bytes(&v)).is_err());
}
#[test]
fn absent_orientation_or_multipart_order_is_a_named_loss_not_a_guessed_feature() {
    let mut v: Value = serde_json::from_slice(HIERARCHY).unwrap();
    object_mut(&mut v, "https://example.org/lab/origin_left")
        .as_object_mut()
        .unwrap()
        .remove("sbol:order");
    let d = import(&bytes(&v)).unwrap();
    let p = d
        .project(&selected(), ProjectionPolicy::AllowOmissions)
        .unwrap();
    assert_eq!(p.molecule.features.len(), 1);
    assert!(p.issues.iter().any(|i| i.message.contains("ordering")));
    let mut v: Value = serde_json::from_slice(HIERARCHY).unwrap();
    object_mut(&mut v, "https://example.org/lab/wrapped")
        .as_object_mut()
        .unwrap()
        .remove("sbol:orientation");
    let d = import(&bytes(&v)).unwrap();
    assert!(
        d.project(&selected(), ProjectionPolicy::AllowOmissions)
            .unwrap()
            .issues
            .iter()
            .any(|i| i.message.contains("unspecified"))
    );
    let mut v: Value = serde_json::from_slice(HIERARCHY).unwrap();
    object_mut(&mut v, "https://example.org/lab/plasmid")["sbol:type"] =
        json!([{"@id":DNA},{"@id":CIRCULAR},{"@id":"http://identifiers.org/SO:0000984"}]);
    assert!(
        import(&bytes(&v))
            .unwrap()
            .project(&selected(), ProjectionPolicy::AllowOmissions)
            .is_err()
    );
}
#[test]
fn standalone_molecule_export_has_stable_iris_and_reverse_join_order_and_explicit_losses() {
    let mut m = Molecule::dna("Selected literal", "AAAACCCCGGGGTTTT", Topology::Circular).unwrap();
    m.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::Complement(Box::new(Location::Join(vec![
            Location::range(0, 4),
            Location::range(8, 12),
        ]))),
        qualifiers: vec![("label".into(), Some("Reverse ordered feature".into()))],
        needs_review: false,
    });
    let ids = ExportIds {
        namespace: "https://example.org/export".into(),
        component: "https://example.org/export/plasmid".into(),
        sequence: "https://example.org/export/sequence".into(),
    };
    let exported = export_molecule(&m, &ids, ProjectionPolicy::Strict).unwrap();
    assert!(exported.issues.is_empty());
    assert_eq!(
        export_molecule(&m, &ids, ProjectionPolicy::Strict)
            .unwrap()
            .bytes,
        exported.bytes
    );
    let doc = import(&exported.bytes).unwrap();
    let converted = doc
        .project(
            &Selection {
                component: ids.component.clone(),
                sequence: ids.sequence.clone(),
                topology: None,
            },
            ProjectionPolicy::Strict,
        )
        .unwrap();
    assert_eq!(
        converted.molecule.features[0]
            .location
            .extract(&converted.molecule.sequence)
            .unwrap(),
        m.features[0].location.extract(&m.sequence).unwrap()
    );
    let original = doc.original_bytes().to_vec();
    m.metadata.insert(
        "duplex_model".into(),
        "foreign chemistry is not transferable authority".into(),
    );
    assert!(export_molecule(&m, &ids, ProjectionPolicy::Strict).is_err());
    let limited = export_molecule(&m, &ids, ProjectionPolicy::AllowOmissions).unwrap();
    assert!(!limited.issues.is_empty());
    assert!(
        !String::from_utf8(limited.bytes)
            .unwrap()
            .contains("foreign chemistry")
    );
    assert_eq!(doc.original_bytes(), original);
}
#[test]
fn bounds_and_missing_sequences_fail_without_network_or_hierarchy_synthesis() {
    assert!(import(&vec![b' '; MAX_BYTES + 1]).is_err());
    let deep = format!("{}0{}", "[".repeat(40), "]".repeat(40));
    assert!(import(deep.as_bytes()).is_err());
    let mut v: Value = serde_json::from_slice(HIERARCHY).unwrap();
    v["@graph"]
        .as_array_mut()
        .unwrap()
        .retain(|o| o["@id"] != "https://example.org/lab/plasmid_sequence");
    let doc = import(&bytes(&v)).unwrap();
    assert!(doc.issues().iter().any(|i| i.code == "unresolved-sequence"));
    assert!(
        doc.project(&selected(), ProjectionPolicy::AllowOmissions)
            .is_err()
    );
    let mut v: Value = serde_json::from_slice(HIERARCHY).unwrap();
    object_mut(&mut v, "https://example.org/lab/plasmid")["sbol:type"] = json!({"@id":DNA});
    let doc = import(&bytes(&v)).unwrap();
    assert!(
        doc.project(&selected(), ProjectionPolicy::AllowOmissions)
            .is_err()
    );
    let mut choice = selected();
    choice.topology = Some(Topology::Linear);
    assert!(
        doc.project(&choice, ProjectionPolicy::AllowOmissions)
            .is_ok()
    );
}

#[test]
fn cut_and_entire_sequence_are_typed_and_retained_without_fabricating_base_ranges() {
    let mut v: Value = serde_json::from_slice(HIERARCHY).unwrap();
    object_mut(&mut v, "https://example.org/lab/wrapped")["sbol:hasLocation"] =
        json!({"@id":"https://example.org/lab/origin_right"});
    v["@graph"]
        .as_array_mut()
        .unwrap()
        .retain(|o| o["@id"] != "https://example.org/lab/origin_left");
    let loc = object_mut(&mut v, "https://example.org/lab/origin_right")
        .as_object_mut()
        .unwrap();
    loc.insert("@type".into(), json!("sbol:Cut"));
    loc.remove("sbol:start");
    loc.remove("sbol:end");
    loc.insert("sbol:at".into(), json!(0));
    let doc = import(&bytes(&v)).unwrap();
    assert!(matches!(
        &doc.objects()["https://example.org/lab/origin_right"].kind,
        Kind::Location(LocationRecord {
            kind: LocationKind::Cut { at: 0 },
            ..
        })
    ));
    let projected = doc
        .project(&selected(), ProjectionPolicy::AllowOmissions)
        .unwrap();
    assert_eq!(projected.molecule.features.len(), 1);
    assert!(
        projected
            .issues
            .iter()
            .any(|i| i.message.contains("Cut boundary"))
    );
    let loc = object_mut(&mut v, "https://example.org/lab/origin_right")
        .as_object_mut()
        .unwrap();
    loc.insert("@type".into(), json!("sbol:EntireSequence"));
    loc.remove("sbol:at");
    let doc = import(&bytes(&v)).unwrap();
    let projected = doc
        .project(&selected(), ProjectionPolicy::AllowOmissions)
        .unwrap();
    assert_eq!(
        projected.molecule.features[1]
            .location
            .extract(&projected.molecule.sequence)
            .unwrap(),
        "AAAACCCCGGGGTTTT"
    );
    let original = doc.original_bytes().to_vec();
    assert_eq!(
        import(&doc.expanded_bytes().unwrap()).unwrap().objects(),
        doc.objects()
    );
    assert_eq!(doc.original_bytes(), original);
}

#[test]
fn explicit_so_type_aliases_preserve_graph_iris_and_semantic_orientation_without_weakening_ssdna_guards()
 {
    const SOURCE: &[u8] = include_bytes!("fixtures/so_aliases.jsonld");
    let selected = Selection {
        component: "https://example.org/aliases/component".into(),
        sequence: "https://example.org/aliases/sequence".into(),
        topology: None,
    };
    let graph = import(SOURCE).unwrap();
    let p = graph.project(&selected, ProjectionPolicy::Strict).unwrap();
    assert_eq!(p.molecule.topology, Topology::Circular);
    assert_eq!(
        p.molecule.features[0]
            .location
            .extract(&p.molecule.sequence)
            .unwrap(),
        "AAAA"
    );
    assert_eq!(
        p.molecule.features[1]
            .location
            .extract(&p.molecule.sequence)
            .unwrap(),
        "GGGG"
    );
    assert_eq!(graph.original_bytes(), SOURCE);
    let expanded = graph.expanded_bytes().unwrap();
    assert!(
        String::from_utf8(expanded.clone())
            .unwrap()
            .contains("https://identifiers.org/SO:0000988")
    );
    assert_eq!(import(&expanded).unwrap().objects(), graph.objects());
    for scheme in ["http", "https"] {
        let mut value: Value = serde_json::from_slice(SOURCE).unwrap();
        object_mut(&mut value, &selected.component)["sbol:type"] = json!([{"@id":DNA},{"@id":format!("{scheme}://identifiers.org/SO:0000987")},{"@id":format!("{scheme}://identifiers.org/SO:0000985")}]);
        assert_eq!(
            import(&bytes(&value))
                .unwrap()
                .project(&selected, ProjectionPolicy::Strict)
                .unwrap()
                .molecule
                .topology,
            Topology::Linear
        );
        object_mut(&mut value, &selected.component)["sbol:type"] = json!([{"@id":DNA},{"@id":CIRCULAR},{"@id":format!("{scheme}://identifiers.org/SO:0000984")}]);
        assert!(
            import(&bytes(&value))
                .unwrap()
                .project(&selected, ProjectionPolicy::AllowOmissions)
                .is_err()
        );
        object_mut(&mut value, &selected.component)["sbol:type"] = json!([{"@id":RNA},{"@id":format!("{scheme}://identifiers.org/SO:0000987")},{"@id":format!("{scheme}://identifiers.org/SO:0000984")}]);
        object_mut(&mut value, &selected.sequence)["sbol:elements"] = json!("ACGUACGU");
        assert_eq!(
            import(&bytes(&value))
                .unwrap()
                .project(&selected, ProjectionPolicy::Strict)
                .unwrap()
                .molecule
                .sequence
                .alphabet(),
            Alphabet::Rna
        );
    }
    let mut conflict: Value = serde_json::from_slice(SOURCE).unwrap();
    object_mut(&mut conflict, &selected.component)["sbol:type"] =
        json!([{"@id":DNA},{"@id":LINEAR},{"@id":"https://identifiers.org/SO:0000988"}]);
    assert!(import(&bytes(&conflict)).is_err());
    object_mut(&mut conflict, &selected.component)["sbol:type"] = json!([{"@id":DNA},{"@id":CIRCULAR},{"@id":"http://identifiers.org/SO:0000984"},{"@id":"https://identifiers.org/SO:0000985"}]);
    assert!(import(&bytes(&conflict)).is_err());
    let mut conflict: Value = serde_json::from_slice(SOURCE).unwrap();
    object_mut(&mut conflict, "https://example.org/aliases/forward_range")["sbol:orientation"] =
        json!({"@id":REVERSE});
    assert!(
        import(&bytes(&conflict))
            .unwrap()
            .project(&selected, ProjectionPolicy::Strict)
            .is_err()
    );
    let plain = Molecule::dna("Canonical export", "AAAACCCC", Topology::Circular).unwrap();
    let written = export_molecule(
        &plain,
        &ExportIds {
            namespace: "https://example.org/out".into(),
            component: "https://example.org/out/component".into(),
            sequence: "https://example.org/out/sequence".into(),
        },
        ProjectionPolicy::Strict,
    )
    .unwrap();
    let text = String::from_utf8(written.bytes).unwrap();
    assert!(text.contains(CIRCULAR));
    assert!(!text.contains("https://identifiers.org/SO:0000988"));
}
