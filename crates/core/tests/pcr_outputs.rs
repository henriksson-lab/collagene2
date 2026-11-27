use collagene_core::{
    analysis::Primer,
    oligo_formats::{self, Delimiter, TableExportMode},
    pcr_outputs::{self, Options, Request},
    protocol::{self, AnalysisJob, AnalysisResult},
    verification::simulate_pcr_documents,
    *,
};
fn fixture() -> Request {
    let mut source = Molecule::dna("template", "ACGACCCCACGAGGGGTAAC", Topology::Linear).unwrap();
    source.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::Complement(Box::new(Location::range(16, 20))),
        qualifiers: vec![("label".into(), Some("reverse site".into()))],
        needs_review: false,
    });
    let forward = Primer::new("=Forward, retained", "ACGA", "CGT").unwrap();
    let reverse = Primer::new("Reverse named", "GTTA", "GAC").unwrap();
    let products = simulate_pcr_documents(&source, &forward, &reverse).unwrap();
    assert_eq!(products.len(), 2);
    assert_eq!(products[1].geometry.forward_start, 8);
    Request {
        source,
        forward,
        reverse,
        selected_index: 1,
        expected_product_revision: products[1].molecule.revision(),
        options: Options {
            name: "Chosen PCR reference".into(),
            save_product: true,
            save_primers: true,
        },
    }
}
#[test]
fn explicit_alternative_keeps_annotation_and_named_complete_oligos() {
    let request = fixture();
    let before = request.source.revision();
    let result = pcr_outputs::prepare(&request).unwrap();
    let product = result.product.as_ref().unwrap();
    assert_eq!(product.name, "Chosen PCR reference");
    assert_eq!(product.sequence.as_str(), "CGTACGAGGGGTAACGTC");
    assert_eq!(
        product.features[0].location,
        Location::Complement(Box::new(Location::range(11, 15)))
    );
    assert_eq!(
        product.features[0]
            .location
            .extract(&product.sequence)
            .unwrap(),
        "GTTA"
    );
    assert_eq!(product.metadata["pcr:template_revision"], before);
    assert!(!product.metadata.contains_key("verification_status"));
    assert!(!product.metadata.contains_key("duplex_model"));
    assert_eq!(request.source.revision(), before);
    let table = result.primers.as_ref().unwrap();
    assert_eq!(table.records[0].name, request.forward.name);
    let f = table.records[0].primer.as_ref().unwrap();
    assert_eq!(f.ordered_sequence(), "CGTACGA");
    assert_eq!(f.annealing.as_str(), "ACGA");
    assert_eq!(f.tail, "CGT");
    assert_eq!(
        table.records[1].primer.as_ref().unwrap().ordered_sequence(),
        "GACGTTA"
    );
    let csv =
        oligo_formats::export_oligos(table, Delimiter::Comma, TableExportMode::SpreadsheetSafe)
            .unwrap();
    assert!(csv.text.contains("'=Forward, retained"));
    assert!(csv.text.contains("CGTACGA"));
    result
        .validate_binding(&pcr_outputs::revision(&request).unwrap())
        .unwrap();
    let mut tampered = result.clone();
    tampered.geometry.sequence = Sequence::dna("ACGT").unwrap();
    assert!(tampered.validate_binding(&result.context_revision).is_err());
    assert!(result.validate_binding(&"0".repeat(64)).is_err());
}
#[test]
fn stale_preview_invalid_geometry_and_unselected_outputs_are_rejected_or_omitted() {
    let request = fixture();
    let mut changed = request.clone();
    changed.selected_index = 0;
    assert!(pcr_outputs::prepare(&changed).is_err());
    changed = request.clone();
    changed.selected_index = 2;
    assert!(pcr_outputs::prepare(&changed).is_err());
    changed = request.clone();
    changed.source.name = "another source snapshot".into();
    assert!(pcr_outputs::prepare(&changed).is_err());
    changed = request.clone();
    changed.forward.name = "different oligo name".into();
    assert!(pcr_outputs::prepare(&changed).is_err());
    changed = request.clone();
    changed.options.save_product = false;
    let only_primers = pcr_outputs::prepare(&changed).unwrap();
    assert!(only_primers.product.is_none());
    assert!(only_primers.primers.is_some());
    changed = request.clone();
    changed.options.save_primers = false;
    let only_product = pcr_outputs::prepare(&changed).unwrap();
    assert!(only_product.product.is_some());
    assert!(only_product.primers.is_none());
    changed.options.save_product = false;
    assert!(pcr_outputs::revision(&changed).is_err());
    changed = request.clone();
    changed.options.name = " invalid ".into();
    assert!(pcr_outputs::revision(&changed).is_err());
    changed = request.clone();
    changed.forward.tail = "Z".into();
    assert!(pcr_outputs::revision(&changed).is_err());
    changed = request.clone();
    changed.expected_product_revision = "bad".into();
    assert!(pcr_outputs::revision(&changed).is_err());
}
fn worker_request(request: Request) -> protocol::Request {
    protocol::Request {
        id: 50,
        revision: pcr_outputs::revision(&request).unwrap(),
        job: AnalysisJob::PreparePcrOutputs { request },
    }
}
#[test]
fn worker_binds_output_selection_and_rejects_forged_old_variant_prefixes() {
    let request = worker_request(fixture());
    let wire = protocol::encode_request(&request).unwrap();
    let mut old = wire.clone();
    old[4..6].copy_from_slice(&49u16.to_le_bytes());
    assert!(matches!(
        protocol::process_request(&old),
        Err(Error::Unsupported(_))
    ));
    // A body truncated after its borrowed prefix must fail the schema gate before
    // decoding any new payload. The unchanged checksum is deliberately irrelevant.
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&wire[42..]).unwrap();
    assert_eq!(tag, 115);
    let payload = postcard::to_allocvec(&(request.id, request.revision.as_str(), 115u32)).unwrap();
    let mut prefix = wire[..42].to_vec();
    prefix[4..6].copy_from_slice(&49u16.to_le_bytes());
    prefix[6..10].copy_from_slice(&(payload.len() as u32).to_le_bytes());
    prefix.extend_from_slice(&payload);
    assert!(matches!(
        protocol::process_request(&prefix),
        Err(Error::Unsupported(_))
    ));
    let response_wire = protocol::process_request(&wire).unwrap();
    let response = protocol::decode_response(&response_wire).unwrap();
    assert_eq!(response.revision, request.revision);
    let AnalysisResult::PcrOutputs(result) = response.result else {
        panic!("wrong response")
    };
    assert_eq!(result.geometry.forward_start, 8);
    let ((_, _, tag), _) =
        postcard::take_from_bytes::<(u32, &str, u32)>(&response_wire[42..]).unwrap();
    assert_eq!(tag, 107);
    let mut old_response = response_wire;
    old_response[4..6].copy_from_slice(&49u16.to_le_bytes());
    assert!(matches!(
        protocol::decode_response(&old_response),
        Err(Error::Unsupported(_))
    ));
    let mut stale = worker_request(fixture());
    stale.revision = "0".repeat(64);
    assert!(protocol::process_request(&protocol::encode_request(&stale).unwrap()).is_err());
}

#[test]
fn iupac_tail_references_and_long_product_only_oligos_keep_existing_model() {
    let mut request = fixture();
    request.forward.tail = "NRYS".into();
    let previews =
        simulate_pcr_documents(&request.source, &request.forward, &request.reverse).unwrap();
    request.expected_product_revision = previews[1].molecule.revision();
    let prepared = pcr_outputs::prepare(&request).unwrap();
    assert_eq!(
        prepared.product.unwrap().sequence.as_str(),
        "NRYSACGAGGGGTAACGTC"
    );
    let table = prepared.primers.unwrap();
    assert_eq!(table.records[0].primer.as_ref().unwrap().tail, "NRYS");
    assert_eq!(
        table.records[0].primer.as_ref().unwrap().ordered_sequence(),
        "NRYSACGA"
    );
    // Long oligos remain order-table evidence even when outside the existing
    // 200-base reusable-primer import profile.
    request.forward.tail = "A".repeat(250);
    let previews =
        simulate_pcr_documents(&request.source, &request.forward, &request.reverse).unwrap();
    request.expected_product_revision = previews[1].molecule.revision();
    let table = pcr_outputs::prepare(&request).unwrap().primers.unwrap();
    assert!(table.records[0].primer.is_none());
    assert_eq!(
        table.records[0].raw_tail.as_deref(),
        Some(request.forward.tail.as_str())
    );
    // A table's existing 8192-byte field limit must not become a PCR limit.
    request.forward.tail = "A".repeat(9000);
    request.options.save_primers = false;
    let previews =
        simulate_pcr_documents(&request.source, &request.forward, &request.reverse).unwrap();
    request.expected_product_revision = previews[1].molecule.revision();
    assert!(
        pcr_outputs::prepare(&request)
            .unwrap()
            .product
            .unwrap()
            .sequence
            .len()
            > 9000
    );
    request.options.save_primers = true;
    assert!(matches!(
        pcr_outputs::prepare(&request),
        Err(Error::TooLarge)
    ));
}
