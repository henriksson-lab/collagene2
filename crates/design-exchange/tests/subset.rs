use collagene_design_exchange::{self as sbol, subset};
use serde_json::{Value, json};
use std::collections::{BTreeMap, BTreeSet};
const C: &str = "https://example.org/selected";
const S: &str = "https://example.org/sequence";
const ALT: &str = "https://example.org/alternative";
fn graph() -> Value {
    json!({"@context":{"sbol":sbol::SBOL},"@graph":[
    {"@id":C,"@type":"sbol:Component","sbol:hasNamespace":{"@id":"https://example.org/"},"sbol:type":{"@id":sbol::DNA},"sbol:hasSequence":[{"@id":S},{"@id":ALT}]},
    {"@id":S,"@type":"sbol:Sequence","sbol:hasNamespace":{"@id":"https://example.org/"},"sbol:encoding":{"@id":sbol::NUCLEIC_ENCODING},"sbol:elements":"ACGT"},
    {"@id":ALT,"@type":"sbol:Sequence","sbol:hasNamespace":{"@id":"https://example.org/"},"sbol:encoding":{"@id":sbol::NUCLEIC_ENCODING},"sbol:elements":"TTTT"},
    {"@id":"https://private.invalid/unrelated","@type":"https://private.invalid/Note","https://private.invalid/value":"SECRET_SENTINEL"}
    ]})
}
fn request() -> subset::Request {
    subset::Request {
        roots: BTreeMap::from([(C.into(), BTreeSet::from([S.into()]))]),
        alternatives: BTreeMap::new(),
    }
}
fn doc(v: &Value) -> sbol::Document {
    sbol::import(&serde_json::to_vec(v).unwrap()).unwrap()
}
#[test]
fn canonical_subset_preserves_exact_selected_nodes_and_never_omitted_private_values() {
    let v = graph();
    let d = doc(&v);
    let original = d.original_bytes().to_vec();
    let e = subset::build(&d, &request()).unwrap();
    assert_eq!(e.retained_objects, BTreeSet::from([C.into(), S.into()]));
    assert_eq!(e.omissions.unrelated_objects, 2);
    assert_eq!(e.omissions.alternative_sequence_edges, 1);
    let text = std::str::from_utf8(&e.bytes).unwrap();
    assert!(!text.contains("SECRET_SENTINEL"));
    assert!(!text.contains(ALT));
    assert!(!text.contains("private.invalid"));
    let received = subset::validate_received(&e.bytes, &request()).unwrap();
    assert_eq!(received.objects().len(), 2);
    assert_eq!(d.original_bytes(), original);
    assert_eq!(subset::build(&d, &request()).unwrap().bytes, e.bytes);
    let mut pretty = e.bytes.clone();
    pretty.push(b' ');
    assert!(subset::validate_received(&pretty, &request()).is_err());
    assert!(subset::validate_received(&original, &request()).is_err());
}
#[test]
fn required_structural_edges_cannot_be_silently_redirected_to_selected_sequence() {
    let mut v = graph();
    v["@graph"][0]["sbol:hasFeature"] = json!({"@id":"https://example.org/feature"});
    v["@graph"].as_array_mut().unwrap().extend([
 json!({"@id":"https://example.org/feature","@type":"sbol:SequenceFeature","sbol:hasLocation":{"@id":"https://example.org/location"}}),
 json!({"@id":"https://example.org/location","@type":"sbol:Range","sbol:hasSequence":{"@id":ALT},"sbol:start":1,"sbol:end":2})]);
    assert!(subset::build(&doc(&v), &request()).is_err());
}
#[test]
fn unsupported_selected_properties_and_classes_and_unused_choices_are_rejected() {
    let mut v = graph();
    v["@graph"][0]["https://private.invalid/reference"] =
        json!({"@id":"https://private.invalid/unrelated"});
    assert!(subset::build(&doc(&v), &request()).is_err());
    let mut v = graph();
    v["@graph"][0]["@type"] = json!(["sbol:Component", "https://example.org/Extra"]);
    assert!(subset::build(&doc(&v), &request()).is_err());
    let mut r = request();
    r.alternatives.insert(
        "https://example.org/not-reached".into(),
        BTreeSet::from([S.into()]),
    );
    assert!(subset::build(&doc(&graph()), &r).is_err());
}
#[test]
fn nested_alternatives_require_explicit_choice_without_first_sequence_default() {
    let mut v = graph();
    let child = "https://example.org/child";
    v["@graph"][0]["sbol:hasFeature"] = json!({"@id":"https://example.org/sub"});
    v["@graph"].as_array_mut().unwrap().extend([
 json!({"@id":"https://example.org/sub","@type":"sbol:SubComponent","sbol:instanceOf":{"@id":child}}),
 json!({"@id":child,"@type":"sbol:Component","sbol:hasNamespace":{"@id":"https://example.org/"},"sbol:type":{"@id":sbol::DNA},"sbol:hasSequence":[{"@id":"https://example.org/childA"},{"@id":"https://example.org/childB"}]}),
 json!({"@id":"https://example.org/childA","@type":"sbol:Sequence","sbol:hasNamespace":{"@id":"https://example.org/"},"sbol:encoding":{"@id":sbol::NUCLEIC_ENCODING},"sbol:elements":"AC"}),
 json!({"@id":"https://example.org/childB","@type":"sbol:Sequence","sbol:hasNamespace":{"@id":"https://example.org/"},"sbol:encoding":{"@id":sbol::NUCLEIC_ENCODING},"sbol:elements":"GT"})]);
    let d = doc(&v);
    assert!(subset::build(&d, &request()).is_err());
    let mut r = request();
    r.alternatives.insert(
        child.into(),
        BTreeSet::from(["https://example.org/childB".into()]),
    );
    let e = subset::build(&d, &r).unwrap();
    assert!(subset::validate_received(&e.bytes, &r).is_ok());
    assert!(!std::str::from_utf8(&e.bytes).unwrap().contains("childA"));
}
#[test]
fn omitted_internal_vocabulary_references_are_not_disclosed() {
    let mut v = graph();
    v["@graph"][0]["sbol:role"] = json!({"@id":"https://private.invalid/unrelated"});
    assert!(subset::build(&doc(&v), &request()).is_err());
}
#[test]
fn independently_selected_root_can_retain_a_shared_alternative_sequence() {
    let mut v = graph();
    let other = "https://example.org/other";
    v["@graph"].as_array_mut().unwrap().push(json!({"@id":other,"@type":"sbol:Component","sbol:hasNamespace":{"@id":"https://example.org/"},"sbol:type":{"@id":sbol::DNA},"sbol:hasSequence":{"@id":ALT}}));
    let mut r = request();
    r.roots.insert(other.into(), BTreeSet::from([ALT.into()]));
    let e = subset::build(&doc(&v), &r).unwrap();
    let got = subset::validate_received(&e.bytes, &r).unwrap();
    assert!(got.objects().contains_key(ALT));
    let sbol::Kind::Component(c) = &got.objects()[C].kind else {
        panic!()
    };
    assert_eq!(c.sequences, vec![S.to_string()]);
}
