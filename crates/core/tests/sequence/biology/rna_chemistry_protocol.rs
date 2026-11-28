use collagene_core::{
    protocol::{
        self,
        request::Request,
        response::Response,
        schema::{AnalysisJob, AnalysisResult},
    },
    rna_chemistry as rna,
};
use collagene_model::{Alphabet, Feature, Location, Molecule, Sequence, Topology};
fn inputs() -> (Molecule, rna::Declaration) {
    let m = Molecule {
        name: "Imported RNA".into(),
        sequence: Sequence::new(Alphabet::Rna, "AC").unwrap(),
        topology: Topology::Linear,
        features: vec![Feature {
            kind: "modified_base".into(),
            location: Location::range(0, 1),
            qualifiers: vec![("mod_base".into(), Some("m6a".into()))],
            needs_review: false,
        }],
        metadata: Default::default(),
        original_genbank: None,
    };
    let candidate = rna::annotation_candidates(&m).unwrap().remove(0);
    let d = rna::Declaration {
        version: rna::DECLARATION_VERSION,
        source_revision: rna::source_revision(&m).unwrap(),
        coverage: rna::Coverage::ClosedModel,
        ends: rna::Ends::Linear {
            five: rna::FivePrime::Hydroxyl,
            three: rna::ThreePrime::Hydroxyl,
        },
        sites: vec![candidate.suggested_site.unwrap()],
        links: vec![],
        annotation_reviews: vec![rna::FeatureReview {
            feature_index: 0,
            expected_feature_revision: candidate.feature_revision,
            site_entry: 0,
        }],
        provenance: rna::Provenance {
            statement: "Explicit fixture hypothesis".into(),
            reference: None,
        },
    };
    (m, d)
}
#[test]
fn binary_and_session_analysis_keep_exact_source_and_declared_chemistry() {
    let (m, d) = inputs();
    let before = m.clone();
    let q = Request {
        id: 31,
        revision: rna::request_revision(&m, &d).unwrap(),
        job: AnalysisJob::AnalyzeRnaChemistry {
            molecule: m.clone(),
            declaration: d.clone(),
        },
    };
    let bytes = protocol::request::Request::encode(&q).unwrap();
    assert_eq!(
        &bytes[4..6],
        &protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes()
    );
    let ((_, _, tag), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(tag, 193);
    let mut session = protocol::session::AnalysisSession::default();
    for response in [
        protocol::execution::process_request(&bytes).unwrap(),
        session.process_request(&bytes).unwrap(),
    ] {
        let ((_, _, tag), _) =
            postcard::take_from_bytes::<(u32, &str, u32)>(&response[42..]).unwrap();
        assert_eq!(tag, 183);
        let decoded = protocol::response::Response::decode(&response).unwrap();
        let AnalysisResult::RnaChemistry(r) = decoded.result else {
            panic!("RNA result")
        };
        r.validate_binding(&q.revision).unwrap();
        r.validate_for(&m, &d).unwrap();
        assert_eq!(r.formula.unwrap().hill(), "C20H27N8O11P");
        assert!(r.formula_complete && r.unreviewed_evidence.is_empty());
        assert_eq!(r.retained_evidence.len(), 1);
    }
    assert_eq!(m, before);
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&115u16.to_le_bytes());
    assert!(protocol::execution::process_request(&old).is_err());
    assert!(session.process_request(&old).is_err());
    let mut stale = q.clone();
    stale.revision = "f".repeat(64);
    assert!(
        protocol::execution::process_request(&protocol::request::Request::encode(&stale).unwrap())
            .is_err()
    );
    let mut invalid = q;
    let AnalysisJob::AnalyzeRnaChemistry { declaration, .. } = &mut invalid.job else {
        unreachable!()
    };
    declaration.sites[0].nucleoside = rna::Nucleoside::Canonical;
    assert!(protocol::request::Request::encode(&invalid).is_err());
}
#[test]
fn checksummed_worker_reports_reject_changed_binding_quantities_and_review_evidence() {
    let (m, d) = inputs();
    let report = rna::analyze(&m, &d).unwrap();
    let refuse = |r: rna::Report| {
        let response = Response {
            id: 31,
            revision: report.context_revision.clone(),
            result: AnalysisResult::RnaChemistry(r),
        };
        let body = postcard::to_allocvec(&response).unwrap();
        let mut bytes = b"CGR\0".to_vec();
        bytes.extend_from_slice(&protocol::schema::WORKER_SCHEMA_VERSION.to_le_bytes());
        bytes.extend_from_slice(&(body.len() as u32).to_le_bytes());
        bytes.extend_from_slice(blake3::hash(&body).as_bytes());
        bytes.extend_from_slice(&body);
        assert!(protocol::response::Response::decode(&bytes).is_err());
    };
    let mut r = report.clone();
    r.mass.as_mut().unwrap().monoisotopic_da += 1.0;
    refuse(r);
    let mut r = report.clone();
    r.formula.as_mut().unwrap().p += 1;
    refuse(r);
    let mut r = report.clone();
    r.unreviewed_evidence = r.retained_evidence.clone();
    refuse(r);
    let mut r = report.clone();
    r.catalog_revision = "f".repeat(64);
    refuse(r);
    let mut r = report.clone();
    r.declaration.annotation_reviews.clear();
    refuse(r);
    let mut d = d;
    d.annotation_reviews.clear();
    d.coverage = rna::Coverage::Partial;
    let r = rna::analyze(&m, &d).unwrap();
    r.validate_binding(&r.context_revision).unwrap();
    assert!(r.mass.is_none());
}
