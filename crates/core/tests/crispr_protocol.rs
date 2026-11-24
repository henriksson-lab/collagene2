//! F25 worker bindings and append-only CRISPR computational-design version gates.
use collagene_core::{
    Error, Molecule, Topology,
    crispr::*,
    edit_template,
    protocol::{self, AnalysisJob, AnalysisResult, Request, Response},
};
const SPACER: &str = "ACGTACGTACGTACGTACGT";
fn guide_fixture() -> (Molecule, GuideOptions) {
    (
        Molecule::dna(
            "reference",
            &format!("{SPACER}AGGTTT{SPACER}AGG"),
            Topology::Linear,
        )
        .unwrap(),
        GuideOptions::default(),
    )
}
fn wire(magic: &[u8; 4], version: u16, value: &impl serde::Serialize) -> Vec<u8> {
    raw(magic, version, &postcard::to_allocvec(value).unwrap())
}
fn raw(magic: &[u8; 4], version: u16, p: &[u8]) -> Vec<u8> {
    let mut b = magic.to_vec();
    b.extend(version.to_le_bytes());
    b.extend((p.len() as u32).to_le_bytes());
    b.extend(blake3::hash(p).as_bytes());
    b.extend(p);
    b
}
fn run(request: &Request) -> Response {
    let bytes = protocol::encode_request(request).unwrap();
    assert_eq!(&bytes[4..6], &23u16.to_le_bytes());
    protocol::decode_response(&protocol::process_request(&bytes).unwrap()).unwrap()
}
#[test]
fn guide_and_selected_reference_worker_reports_bind_exact_roles_labels_options_and_coordinates() {
    let (molecule, options) = guide_fixture();
    let request = Request {
        id: 601,
        revision: guide_revision(&molecule, &options).unwrap(),
        job: AnalysisJob::DesignCrisprGuides {
            molecule: molecule.clone(),
            options: options.clone(),
        },
    };
    let AnalysisResult::CrisprGuides(report) = run(&request).result else {
        panic!()
    };
    report.validate_bindings_for(&molecule, &options).unwrap();
    let candidate = report
        .candidates
        .iter()
        .find(|g| !g.reverse && g.guide_positions[0] == 0 && g.guide_sequence == SPACER)
        .unwrap();
    validate_candidate(&molecule, &options, candidate).unwrap();
    let refs = vec![
        ReferenceInput {
            id: "intended".into(),
            molecule: molecule.clone(),
            build: Some("declared-build".into()),
            accession: Some("ref.1".into()),
        },
        ReferenceInput {
            id: "same-sequence-other-role".into(),
            molecule: molecule.clone(),
            build: None,
            accession: None,
        },
    ];
    let query = GuideQuery {
        id: candidate.id.clone(),
        sequence: candidate.guide_sequence.clone(),
    };
    let mut off = OffTargetOptions::default();
    off.max_mismatches = 0;
    off.intended = Some(IntendedTarget {
        reference_id: "intended".into(),
        reference_revision: molecule.revision(),
        guide_positions: candidate.guide_positions.clone(),
        pam_positions: candidate.pam_positions.clone(),
        reverse: candidate.reverse,
    });
    let revision = off_target_revision(&query, &refs, &off).unwrap();
    let request = Request {
        id: 602,
        revision: revision.clone(),
        job: AnalysisJob::SearchCrisprOffTargets {
            query: query.clone(),
            references: refs.clone(),
            options: off.clone(),
        },
    };
    let AnalysisResult::CrisprOffTargets(report) = run(&request).result else {
        panic!()
    };
    report.validate_bindings_for(&query, &refs, &off).unwrap();
    replay_off_targets(&query, &refs, &off, &report).unwrap();
    assert_eq!(report.hits.iter().filter(|h| h.is_intended).count(), 1);
    assert_eq!(report.hits.len(), 4);
    let mut changed = refs.clone();
    changed[0].build = Some("another-label".into());
    assert_ne!(
        off_target_revision(&query, &changed, &off).unwrap(),
        revision
    );
    assert!(
        report
            .validate_bindings_for(&query, &changed, &off)
            .is_err()
    );
    changed = refs.clone();
    changed.swap(0, 1);
    assert!(
        report
            .validate_bindings_for(&query, &changed, &off)
            .is_err()
    );
    let mut stale = request;
    stale.revision = "0".repeat(64);
    assert!(protocol::process_request(&protocol::encode_request(&stale).unwrap()).is_err());
    let response = Response {
        id: 602,
        revision: "0".repeat(64),
        result: AnalysisResult::CrisprOffTargets(report),
    };
    assert!(protocol::decode_response(&wire(b"CGR\0", 18, &response)).is_err());
}
fn donor_fixture() -> (Molecule, Molecule, edit_template::EditTemplateOptions) {
    (
        Molecule::dna("source", "AAAACCCCGGGGTTTT", Topology::Linear).unwrap(),
        Molecule::dna("target", "AAAACCTCGGGGTTTT", Topology::Linear).unwrap(),
        edit_template::EditTemplateOptions {
            edit_start: 6,
            delete_length: 1,
            replacement: "T".into(),
            left_arm: 3,
            right_arm: 3,
            reverse_donor: false,
            max_cut_distance: 100,
        },
    )
}
#[test]
fn edit_template_design_and_independent_check_bind_supplied_donor_without_physical_edit_authority()
{
    let (source, target, options) = donor_fixture();
    let revision = edit_template::batch_revision(&source, &target, &options, None).unwrap();
    let request = Request {
        id: 603,
        revision,
        job: AnalysisJob::DesignEditTemplate {
            source: source.clone(),
            target: target.clone(),
            options: options.clone(),
            guide: None,
        },
    };
    let AnalysisResult::EditTemplate(report) = run(&request).result else {
        panic!()
    };
    report
        .validate_bindings_for(&source, &target, &options, None)
        .unwrap();
    assert!(report.check.consistent);
    assert_eq!(report.donor.sequence.as_str(), "ACCTCGG");
    assert_eq!(report.edited_reference.sequence, target.sequence);
    let donor = report.donor.sequence.as_str().to_owned();
    let revision =
        edit_template::check_batch_revision(&source, &target, &options, None, &donor).unwrap();
    let request = Request {
        id: 604,
        revision,
        job: AnalysisJob::CheckEditTemplate {
            source: source.clone(),
            target: target.clone(),
            options: options.clone(),
            guide: None,
            donor: donor.clone(),
        },
    };
    let AnalysisResult::EditTemplateCheck(check) = run(&request).result else {
        panic!()
    };
    check
        .validate_bindings_for(&source, &target, &options, None, &donor)
        .unwrap();
    assert!(check.consistent);
    assert!(
        check
            .validate_bindings_for(&source, &target, &options, None, "ACCTCGT")
            .is_err()
    );
    let wrong = edit_template::check(&source, &target, &options, None, "ACCTCGT").unwrap();
    assert!(!wrong.consistent);
    assert!(!wrong.right_matches);
}
#[test]
fn all_four_new_jobs_and_results_require_schema_eighteen_before_body_decode_and_stay_ordinary_sized()
 {
    for (magic, tags) in [
        (b"CGQ\0", vec![49u32, 50, 51, 52]),
        (b"CGR\0", vec![47u32, 48, 49, 50]),
    ] {
        for tag in tags {
            let prefix = postcard::to_allocvec(&(605u32, "batch", tag)).unwrap();
            for version in 1..18 {
                let b = raw(magic, version, &prefix);
                let error = if magic == b"CGQ\0" {
                    protocol::process_request(&b).unwrap_err()
                } else {
                    protocol::decode_response(&b).unwrap_err()
                };
                assert!(error.to_string().contains("schema 18"), "{error}");
            }
            let mut large = prefix;
            large.resize(collagene_core::formats::MAX_IMPORT_BYTES + 1, 0);
            let b = raw(magic, 18, &large);
            let error = if magic == b"CGQ\0" {
                protocol::process_request(&b).unwrap_err()
            } else {
                protocol::decode_response(&b).unwrap_err()
            };
            assert!(matches!(error, Error::TooLarge));
        }
    }
}
#[test]
fn scientific_source_bounds_precede_revision_hashes_and_stale_comparisons() {
    let (mut molecule, options) = guide_fixture();
    molecule
        .metadata
        .insert("oversized evidence".into(), "x".repeat(17 * 1024 * 1024));
    let request = Request {
        id: 606,
        revision: "stale".into(),
        job: AnalysisJob::DesignCrisprGuides { molecule, options },
    };
    assert!(matches!(
        protocol::process_request(&protocol::encode_request(&request).unwrap()),
        Err(Error::TooLarge)
    ));
    let (mut source, target, options) = donor_fixture();
    source
        .metadata
        .insert("oversized evidence".into(), "x".repeat(9 * 1024 * 1024));
    let request = Request {
        id: 607,
        revision: "stale".into(),
        job: AnalysisJob::DesignEditTemplate {
            source,
            target,
            options,
            guide: None,
        },
    };
    assert!(matches!(
        protocol::process_request(&protocol::encode_request(&request).unwrap()),
        Err(Error::TooLarge)
    ));
}
