use collagene_core::{
    cloning::pcr::{
        mismatch_outputs::{self as mismatch_pcr_outputs, Request},
        mutagenesis::PcrOptions,
        outputs::Options,
        primers::BindingOptions,
    },
    formats::oligos::{self as oligo_formats, Delimiter, TableExportMode},
    protocol::{
        self,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::analysis::Primer,
    verification::pcr::simulate_mutagenesis_preview,
};
use collagene_model::{Error, Feature, Location, Molecule, Sequence, Topology};
fn fixture() -> Request {
    let mut source = Molecule::dna(
        "original template",
        "ACGACCCCACGAGGGGTAAC",
        Topology::Linear,
    )
    .unwrap();
    source.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::Complement(Box::new(Location::range(16, 20))),
        qualifiers: vec![("label".into(), Some("changed reverse-binding site".into()))],
        needs_review: false,
    });
    let forward = Primer::new("=Forward, named", "ATGA", "CGT").unwrap();
    let reverse = Primer::new("Reverse named", "GTCA", "GAC").unwrap();
    let options = PcrOptions {
        binding: BindingOptions {
            max_mismatches: 1,
            min_three_prime_exact: 1,
            ..Default::default()
        },
        ..Default::default()
    };
    let previews = simulate_mutagenesis_preview(&source, &forward, &reverse, &options).unwrap();
    assert!(previews.documents.len() >= 2);
    let selected_index = previews
        .report
        .products
        .iter()
        .position(|g| {
            g.forward_start == 8
                && g.reverse_start == 16
                && g.forward_primer_index == 0
                && g.reverse_primer_index == 1
        })
        .unwrap();
    Request {
        source,
        forward,
        reverse,
        options,
        selected_index,
        expected_product_revision: previews.documents[selected_index].revision(),
        outputs: Options {
            name: "Chosen mismatch-PCR reference".into(),
            save_product: true,
            save_primers: true,
        },
    }
}
#[test]
fn selected_alternative_retains_literal_substitutions_annotations_and_named_oligos() {
    let request = fixture();
    let before = request.source.revision();
    let prepared = mismatch_pcr_outputs::prepare(&request).unwrap();
    let product = prepared.product.as_ref().unwrap();
    assert_eq!(product.sequence.as_str(), "CGTATGAGGGGTGACGTC");
    assert_eq!(product.name, "Chosen mismatch-PCR reference");
    assert_eq!(
        product.features[0].location,
        Location::Complement(Box::new(Location::range(11, 15)))
    );
    assert!(product.features[0].needs_review);
    assert_eq!(
        prepared
            .geometry
            .substitutions
            .iter()
            .map(|s| (
                s.source_position,
                s.product_position,
                s.template_base,
                s.product_base
            ))
            .collect::<Vec<_>>(),
        vec![(9, 4, 'C', 'T'), (17, 12, 'A', 'G')]
    );
    assert_eq!(prepared.geometry.forward_binding.mismatches, vec![1]);
    assert_eq!(prepared.geometry.reverse_binding.mismatches, vec![2]);
    assert_eq!(product.metadata["pcr:template_revision"], before);
    assert!(product.metadata["pcr:model"].contains("no target verification"));
    assert!(!product.metadata.contains_key("duplex_model"));
    assert_eq!(request.source.revision(), before);
    let table = prepared.primers.as_ref().unwrap();
    assert_eq!(table.records[0].name, request.forward.name);
    assert_eq!(
        table.records[0].primer.as_ref().unwrap().ordered_sequence(),
        "CGTATGA"
    );
    assert_eq!(
        table.records[1].primer.as_ref().unwrap().ordered_sequence(),
        "GACGTCA"
    );
    let csv =
        oligo_formats::export_oligos(table, Delimiter::Comma, TableExportMode::SpreadsheetSafe)
            .unwrap();
    assert!(csv.text.contains("'=Forward, named"));
    prepared
        .validate_binding(&mismatch_pcr_outputs::revision(&request).unwrap())
        .unwrap();
    let mut tampered = prepared.clone();
    tampered.geometry.sequence = Sequence::dna("ACGT").unwrap();
    assert!(
        tampered
            .validate_binding(&prepared.context_revision)
            .is_err()
    );
}
#[test]
fn stale_selected_preview_and_unselected_output_choices_remain_explicit() {
    let request = fixture();
    let mut changed = request.clone();
    changed.source.name = "another template snapshot".into();
    assert!(mismatch_pcr_outputs::prepare(&changed).is_err());
    changed = request.clone();
    changed.forward.name = "another oligo name".into();
    assert!(mismatch_pcr_outputs::prepare(&changed).is_err());
    changed = request.clone();
    changed.options.binding.min_three_prime_exact = 2;
    assert!(mismatch_pcr_outputs::prepare(&changed).is_err());
    changed = request.clone();
    changed.expected_product_revision = "0".repeat(64);
    assert!(mismatch_pcr_outputs::prepare(&changed).is_err());
    changed = request.clone();
    changed.selected_index = changed.options.max_products;
    assert!(mismatch_pcr_outputs::revision(&changed).is_err());
    changed = request.clone();
    changed.outputs.save_product = false;
    let result = mismatch_pcr_outputs::prepare(&changed).unwrap();
    assert!(result.product.is_none());
    assert!(result.primers.is_some());
    changed = request.clone();
    changed.outputs.save_primers = false;
    let result = mismatch_pcr_outputs::prepare(&changed).unwrap();
    assert!(result.product.is_some());
    assert!(result.primers.is_none());
    changed.outputs.save_product = false;
    assert!(mismatch_pcr_outputs::revision(&changed).is_err());
}
#[test]
fn selected_capped_overlap_preview_does_not_become_a_complete_reaction_claim() {
    let source = Molecule::dna("repeated circle", "ATATATATAT", Topology::Circular).unwrap();
    let forward = Primer::new("one oligo", "ATAT", "").unwrap();
    let reverse = forward.clone();
    let options = PcrOptions {
        binding: BindingOptions::default(),
        max_products: 1,
        ..Default::default()
    };
    let preview = simulate_mutagenesis_preview(&source, &forward, &reverse, &options).unwrap();
    assert_eq!(preview.documents.len(), 1);
    assert!(!preview.report.search_complete);
    assert!(!preview.report.enumeration_complete);
    let request = Request {
        source,
        forward,
        reverse,
        options,
        selected_index: 0,
        expected_product_revision: preview.documents[0].revision(),
        outputs: Options {
            name: "Conditional overlap preview".into(),
            save_product: true,
            save_primers: false,
        },
    };
    let result = mismatch_pcr_outputs::prepare(&request).unwrap();
    assert_eq!(result.geometry.sequence.as_str(), "ATAT");
    assert!(result.geometry.overlapping_domains);
    assert!(!result.search.search_complete);
    assert!(!result.search.enumeration_complete);
    assert!(result.search.unsupported_geometries > 0);
    assert_eq!(result.search.options.max_products, 1);
    assert_eq!(result.search.products_reported, 1);
    assert!(
        result
            .search
            .warnings
            .iter()
            .any(|w| w.contains("truncated"))
    );
    let mut forged = result.clone();
    forged.search.enumeration_complete = true;
    assert!(forged.validate_binding(&result.context_revision).is_err());
}
#[test]
fn worker_51_binds_selection_and_rejects_borrowed_new_tags_under_old_versions() {
    let payload = fixture();
    let revision = mismatch_pcr_outputs::revision(&payload).unwrap();
    let request = protocol::request::Request {
        id: 51,
        revision: revision.clone(),
        job: AnalysisJob::PrepareMismatchPcrOutputs { request: payload },
    };
    let wire = protocol::request::Request::encode(&request).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&wire[42..]).unwrap();
    assert_eq!(tag, 116);
    let mut old = wire.clone();
    old[4..6].copy_from_slice(&50u16.to_le_bytes());
    assert!(matches!(
        protocol::execution::process_request(&old),
        Err(Error::Unsupported(_))
    ));
    let prefix = postcard::to_allocvec(&(request.id, request.revision.as_str(), 116u32)).unwrap();
    let mut truncated = wire[..42].to_vec();
    truncated[4..6].copy_from_slice(&50u16.to_le_bytes());
    truncated[6..10].copy_from_slice(&(prefix.len() as u32).to_le_bytes());
    truncated.extend(prefix);
    assert!(matches!(
        protocol::execution::process_request(&truncated),
        Err(Error::Unsupported(_))
    ));
    let response_wire = protocol::execution::process_request(&wire).unwrap();
    let ((_, _, tag), _) =
        postcard::take_from_bytes::<(u32, &str, u32)>(&response_wire[42..]).unwrap();
    assert_eq!(tag, 108);
    let response = protocol::response::Response::decode(&response_wire).unwrap();
    let AnalysisResult::MismatchPcrOutputs(prepared) = response.result else {
        panic!("wrong output")
    };
    prepared.validate_binding(&revision).unwrap();
    assert_eq!(prepared.geometry.forward_start, 8);
    let mut old_response = response_wire;
    old_response[4..6].copy_from_slice(&50u16.to_le_bytes());
    assert!(matches!(
        protocol::response::Response::decode(&old_response),
        Err(Error::Unsupported(_))
    ));
    let mut stale = request;
    stale.revision = "0".repeat(64);
    assert!(
        protocol::execution::process_request(&protocol::request::Request::encode(&stale).unwrap())
            .is_err()
    );
}
