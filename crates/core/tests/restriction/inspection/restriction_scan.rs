//! Hand-specified recognition evidence; no product reconstruction is a scan oracle.
use collagene_core::{
    analysis::{self, Enzyme},
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    restriction_inspection::*,
};
use collagene_model::{Molecule, Topology};

fn enzyme(name: &str) -> Enzyme {
    analysis::starter_enzymes()
        .into_iter()
        .find(|e| e.name == name)
        .unwrap()
}
fn dna(sequence: &str, topology: Topology) -> Molecule {
    Molecule::dna("literal scan", sequence, topology).unwrap()
}
fn scan_only(report: &RestrictionInspection) {
    assert!(!report.products_determined);
    assert!(report.fragments.is_empty());
    assert!(report.ref_top_lengths.is_none());
    report.validate().unwrap();
}

#[test]
fn literal_ecori_and_outside_bsai_recognition_match_inspection_without_products() {
    let molecule = dna("AAAAGAATTCTTTTGGTCTC", Topology::Linear);
    let enzymes = [enzyme("EcoRI"), enzyme("BsaI")];
    let options = InspectionOptions::default();
    let scanned = scan("source", &molecule, &enzymes, &options).unwrap();
    scan_only(&scanned);
    assert!(scanned.scan_complete);
    assert_eq!(scanned.sites.len(), 2);
    assert_eq!(
        (
            scanned.sites[0].recognition_start,
            scanned.sites[0].top_cut,
            scanned.sites[0].bottom_cut
        ),
        (4, Some(5), Some(9))
    );
    assert!(scanned.sites[0].certain);
    let outside = &scanned.sites[1];
    assert_eq!(
        (
            outside.recognition_start,
            outside.top_boundary,
            outside.bottom_boundary
        ),
        (14, 21, 25)
    );
    assert_eq!((outside.top_cut, outside.bottom_cut), (None, None));
    assert!(
        outside
            .diagnostics
            .contains(&SiteDiagnostic::CutsOutsideSubstrate)
    );
    let inspected = inspect("source", &molecule, &enzymes, &options).unwrap();
    assert_eq!(scanned.sites, inspected.sites);
    assert_eq!(scanned.comparisons, inspected.comparisons);
    assert_eq!(scanned.batch_revision, inspected.batch_revision);
    scanned
        .validate_bindings_for("source", &molecule, &enzymes, &options)
        .unwrap();
    // Even a physically reconstructible case must stay scan-only through this API.
    let exact = dna("TTTGAATTCAAA", Topology::Linear);
    let scanned = scan("source", &exact, &enzymes[..1], &options).unwrap();
    scan_only(&scanned);
    assert!(
        inspect("source", &exact, &enzymes[..1], &options)
            .unwrap()
            .products_determined
    );
}

#[test]
fn origin_wrapping_and_ambiguous_source_bases_keep_reference_evidence_honest() {
    let options = InspectionOptions::default();
    let wrapped = scan(
        "circle",
        &dna("ATTCAAAAGA", Topology::Circular),
        &[enzyme("EcoRI")],
        &options,
    )
    .unwrap();
    scan_only(&wrapped);
    assert!(wrapped.scan_complete);
    assert_eq!(wrapped.sites.len(), 1);
    let site = &wrapped.sites[0];
    assert_eq!(
        (site.recognition_start, site.top_cut, site.bottom_cut),
        (8, Some(9), Some(3))
    );
    assert_eq!(
        site.spans,
        vec![
            ReferenceSpan { start: 8, end: 10 },
            ReferenceSpan { start: 0, end: 4 }
        ]
    );
    let ambiguous = scan(
        "possible",
        &dna("AAAAGNATTCTTTT", Topology::Linear),
        &[enzyme("EcoRI")],
        &options,
    )
    .unwrap();
    scan_only(&ambiguous);
    assert!(ambiguous.scan_complete);
    assert_eq!(ambiguous.sites.len(), 1);
    assert!(!ambiguous.sites[0].certain);
    assert!(
        ambiguous.sites[0]
            .diagnostics
            .contains(&SiteDiagnostic::PossibleRecognition)
    );
    let degenerate = Enzyme {
        name: "literal IUPAC".into(),
        recognition: "GGNCC".into(),
        top_cut: 1,
        bottom_cut: 4,
    };
    let concrete = scan(
        "iupac",
        &dna("GGACC", Topology::Linear),
        &[degenerate],
        &options,
    )
    .unwrap();
    assert_eq!(concrete.sites.len(), 1);
    assert!(concrete.sites[0].certain);
    scan_only(&concrete);
}

#[test]
fn comparison_bound_preserves_observed_site_without_claiming_complete_counts() {
    let molecule = dna("GAATTCGAATTC", Topology::Linear);
    let enzymes = [enzyme("EcoRI")];
    let options = InspectionOptions {
        max_comparisons: 8,
        ..InspectionOptions::default()
    };
    let scanned = scan("bounded", &molecule, &enzymes, &options).unwrap();
    scan_only(&scanned);
    assert_eq!(scanned.comparisons, 8);
    assert!(!scanned.scan_complete);
    assert_eq!(scanned.sites.len(), 1);
    assert_eq!(scanned.sites[0].recognition_start, 0);
    assert!(scanned.limits.contains(&InspectionLimit::Comparisons));
    let inspected = inspect("bounded", &molecule, &enzymes, &options).unwrap();
    assert_eq!(scanned.sites, inspected.sites);
    assert_eq!(scanned.limits, inspected.limits);
}

#[test]
fn current_worker_scan_preserves_scientific_evidence_and_rejects_unsupported_headers() {
    let molecule = dna("TTTGAATTCAAA", Topology::Linear);
    let enzymes = vec![enzyme("EcoRI")];
    let options = InspectionOptions::default();
    let revision = batch_revision("source", &molecule, &enzymes, &options).unwrap();
    let request = Request {
        id: 7,
        revision: revision.clone(),
        job: AnalysisJob::ScanRestriction {
            source_role: "source".into(),
            molecule: molecule.clone(),
            enzymes: enzymes.clone(),
            options: options.clone(),
        },
    };
    let encoded = protocol::request::Request::encode(&request).unwrap();
    assert_eq!(
        &encoded[4..6],
        &collagene_core::protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    let ((id, binding, tag), _) =
        postcard::take_from_bytes::<(u32, &str, u32)>(&encoded[42..]).unwrap();
    assert_eq!((id, binding, tag), (7, revision.as_str(), 76));
    let mut old_request = encoded.clone();
    old_request[4..6].copy_from_slice(&22u16.to_le_bytes());
    assert!(matches!(
        protocol::execution::process_request(&old_request),
        Err(collagene_model::Error::Unsupported(_))
    ));
    let response = protocol::execution::process_request(&encoded).unwrap();
    let decoded = protocol::response::Response::decode(&response).unwrap();
    assert_eq!(decoded.revision, revision);
    let AnalysisResult::RestrictionInspection(report) = decoded.result else {
        panic!("expected scan evidence")
    };
    scan_only(&report);
    assert_eq!(
        &response[4..6],
        &protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    assert!(protocol::response::Response::decode(&response).is_ok());
    let old = Request {
        id: 7,
        revision,
        job: AnalysisJob::InspectRestriction {
            source_role: "source".into(),
            molecule,
            enzymes,
            options,
        },
    };
    let old = protocol::request::Request::encode(&old).unwrap();
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&old[42..]).unwrap();
    assert_eq!(tag, 67);
}
