use collagene_core::{
    analysis::Primer,
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    verification::simulate_pcr_documents,
    *,
};

fn fixture() -> (Molecule, Primer, Primer) {
    let mut source = Molecule::dna("template", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
    source.features.push(Feature {
        kind: "primer_bind".into(),
        location: Location::Complement(Box::new(Location::range(12, 16))),
        qualifiers: vec![("label".into(), Some("reverse binding".into()))],
        needs_review: false,
    });
    (
        source,
        Primer::new("F", "ACGA", "CGT").unwrap(),
        Primer::new("R", "GTTA", "GAC").unwrap(),
    )
}

#[test]
fn standalone_preview_keeps_tail_coordinates_and_immutable_source_binding() {
    let (source, f, r) = fixture();
    let products =
        simulate_pcr_documents(&source, &f, &r, collagene_core::chemistry::Strand::Top).unwrap();
    assert_eq!(products.len(), 1);
    let product = &products[0].molecule;
    assert_eq!(product.sequence.as_str(), "CGTACGACCCCGGGGTAACGTC");
    assert_eq!(
        product.features[0].location,
        Location::Complement(Box::new(Location::range(15, 19)))
    );
    assert_eq!(
        product.features[0]
            .location
            .extract(&product.sequence)
            .unwrap(),
        "GTTA"
    );
    assert_eq!(product.metadata["pcr:template_revision"], source.revision());
    assert!(!product.metadata.contains_key("duplex_model"));
    assert!(!product.metadata.contains_key("verification_status"));
    let saved = protocol::save_document(product).unwrap();
    let reopened = protocol::load_document(&saved).unwrap();
    assert_eq!(reopened.revision(), product.revision());
    assert_eq!(
        source.features[0].location,
        Location::Complement(Box::new(Location::range(12, 16)))
    );
}

#[test]
fn worker_returns_projected_document_instead_of_flattening_sequence() {
    let (source, f, r) = fixture();
    let request = Request {
        id: 127,
        revision: source.revision(),
        job: AnalysisJob::PcrDocuments {
            molecule: source.clone(),
            forward: f,
            reverse: r,
            template_strand: collagene_core::chemistry::Strand::Top,
        },
    };
    let reply = protocol::process_request(&protocol::encode_request(&request).unwrap()).unwrap();
    let response = protocol::decode_response(&reply).unwrap();
    assert_eq!(response.id, 127);
    assert_eq!(response.revision, source.revision());
    let AnalysisResult::PcrDocuments(products) = response.result else {
        panic!("wrong response")
    };
    assert_eq!(
        products[0].molecule.features[0]
            .location
            .extract(&products[0].molecule.sequence)
            .unwrap(),
        "GTTA"
    );
}

#[test]
fn physical_or_declared_single_strand_templates_are_not_flattened() {
    let (mut source, f, r) = fixture();
    source
        .metadata
        .insert("genbank.molecule_type".into(), "ss-DNA".into());
    assert!(matches!(
        simulate_pcr_documents(&source, &f, &r, collagene_core::chemistry::Strand::Top),
        Err(Error::Unsupported(_))
    ));
    source.metadata.remove("genbank.molecule_type");
    source
        .metadata
        .insert("duplex_model".into(), "untrusted attachment".into());
    assert!(matches!(
        simulate_pcr_documents(&source, &f, &r, collagene_core::chemistry::Strand::Top),
        Err(Error::Unsupported(_))
    ));
}

#[test]
fn candidate_schema_one_cannot_claim_new_physical_operations() {
    use collagene_core::{
        transfer::{self, Candidate},
        verification::{DesignPlan, Operation, Source},
    };
    let target = Molecule::dna("target", "ACGT", Topology::Linear).unwrap();
    let candidate = Candidate::new(DesignPlan {
        name: "new physical step".into(),
        target_revision: target.revision(),
        inputs: vec![],
        steps: vec![Operation::Digest {
            template: Source::Input("missing".into()),
            enzymes: vec!["EcoRI".into()],
            fragment_index: 0,
        }],
    });
    let mut bytes = transfer::save_candidate(&candidate).unwrap();
    assert!(transfer::load_candidate(&bytes).is_ok());
    bytes[4..6].copy_from_slice(&1u16.to_le_bytes());
    assert!(matches!(
        transfer::load_candidate(&bytes),
        Err(Error::Unsupported(_))
    ));
}
