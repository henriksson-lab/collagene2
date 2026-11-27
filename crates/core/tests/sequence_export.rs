use collagene_core::{
    Alphabet, Feature, Location, Molecule, Sequence, Topology, gcg_sequence, legacy_formats,
    protocol,
    sequence_export::{self as export, Format},
    snapgene,
};
fn molecule() -> Molecule {
    Molecule::dna("Current reference", "ACGTTAGC", Topology::Linear).unwrap()
}
fn gcg(id: &str, date: &str) -> Format {
    Format::Gcg {
        options: gcg_sequence::ExportOptions {
            record_id: Some(id.into()),
            date_label: Some(date.into()),
        },
    }
}
fn same<T: serde::Serialize>(a: &T, b: &T) {
    assert_eq!(
        serde_json::to_value(a).unwrap(),
        serde_json::to_value(b).unwrap()
    )
}
#[test]
fn snapgene_review_and_prepared_bytes_share_exact_source_and_context_without_spurious_ack() {
    let m = molecule();
    let format = Format::SnapGene;
    let review = export::review(&m, &format).unwrap();
    assert_eq!(review.source_revision, m.revision());
    assert_eq!(
        review.context_revision,
        export::revision(&m, &format).unwrap()
    );
    assert!(review.losses.is_empty());
    assert!(
        !review.warnings.is_empty(),
        "Advisory scope warnings do not invent a lost feature"
    );
    let prepared = export::prepare(&m, &format, &review.context_revision, false).unwrap();
    same(&review, &prepared.assessment);
    prepared.validate_binding(&review.context_revision).unwrap();
    let read = snapgene::import_snapgene(&prepared.bytes)
        .unwrap()
        .molecules
        .remove(0);
    assert_eq!(read.sequence, m.sequence);
    assert_eq!(read.name, m.name);
}
#[test]
fn xdna_name_loss_requires_explicit_current_ack_and_preserves_current_sequence() {
    let m = molecule();
    let review = export::review(&m, &Format::Xdna).unwrap();
    assert!(review.losses.iter().any(|l| l.path == "name"));
    assert!(export::prepare(&m, &Format::Xdna, &review.context_revision, false).is_err());
    let prepared = export::prepare(&m, &Format::Xdna, &review.context_revision, true).unwrap();
    same(&review, &prepared.assessment);
    let read = legacy_formats::import_xdna(&prepared.bytes, &m.name)
        .unwrap()
        .molecules
        .remove(0);
    assert_eq!(read.sequence, m.sequence);
    let mut changed = m.clone();
    changed.name = "Changed exact source identity".into();
    assert!(export::prepare(&changed, &Format::Xdna, &review.context_revision, true).is_err());
    assert!(export::prepare(&m, &Format::SnapGene, &review.context_revision, true).is_err());
}
#[test]
fn gcg_omission_warnings_stay_acknowledged_and_explicit_options_bind_the_context() {
    let mut m = molecule();
    m.topology = Topology::Circular;
    m.metadata
        .insert("description".into(), "Native source detail".into());
    m.features.push(Feature {
        kind: "gene".into(),
        location: Location::range(0, 4),
        qualifiers: vec![("label".into(), Some("Native annotation".into()))],
        needs_review: false,
    });
    let format = gcg("declared_id", "26-NOV-2026");
    let review = export::review(&m, &format).unwrap();
    assert!(review.losses.iter().any(|l| l.message.contains("features")));
    assert!(
        review
            .losses
            .iter()
            .any(|l| l.message.contains("Circular topology"))
    );
    assert!(review.losses.iter().any(|l| l.path == "name"));
    assert!(export::prepare(&m, &format, &review.context_revision, false).is_err());
    let prepared = export::prepare(&m, &format, &review.context_revision, true).unwrap();
    same(&review, &prepared.assessment);
    let read = gcg_sequence::import_gcg(
        &prepared.bytes,
        &gcg_sequence::ImportOptions {
            alphabet: Alphabet::Dna,
            fallback_name: "unused".into(),
        },
    )
    .unwrap()
    .molecules
    .remove(0);
    assert_eq!(read.sequence, m.sequence);
    assert!(read.features.is_empty());
    assert!(
        export::prepare(
            &m,
            &gcg("different_id", "26-NOV-2026"),
            &review.context_revision,
            true
        )
        .is_err()
    );
    assert!(
        export::prepare(
            &m,
            &gcg("declared_id", "27-NOV-2026"),
            &review.context_revision,
            true
        )
        .is_err()
    );
    for format in [
        gcg("bad id", "26-NOV-2026"),
        gcg("id", "Length: 8"),
        Format::Gcg {
            options: Default::default(),
        },
    ] {
        assert!(export::revision(&m, &format).is_err())
    }
}
#[test]
fn physical_and_modified_protein_guards_are_never_overridden_by_acknowledgement() {
    for key in [
        "duplex_model",
        "tracked_modifications",
        "source:declared_termini",
        "source:recombination:reference_only",
    ] {
        let mut m = molecule();
        m.metadata.insert(key.into(), "true".into());
        for format in [Format::SnapGene, Format::Xdna, gcg("id", "26-NOV-2026")] {
            assert!(export::revision(&m, &format).is_err(), "{key}");
            assert!(export::review(&m, &format).is_err());
            assert!(export::prepare(&m, &format, &"0".repeat(64), true).is_err())
        }
    }
    let mut m = molecule();
    m.sequence = Sequence::new(Alphabet::Protein, "MKWVTF").unwrap();
    m.metadata.insert(
        "protein:modifications".into(),
        "retain exact annotation".into(),
    );
    for format in [Format::SnapGene, Format::Xdna, gcg("id", "26-NOV-2026")] {
        assert!(export::review(&m, &format).is_err())
    }
}
#[test]
fn documentary_reply_validation_rejects_bad_context_profile_paths_and_container_shapes() {
    let mut m = molecule();
    m.name = "../CON/unsafe:name".into();
    let review = export::review(&m, &Format::SnapGene).unwrap();
    assert!(!review.filename.contains('/'));
    assert!(!review.filename.contains(':'));
    assert!(!review.filename.starts_with('.'));
    assert!(
        review
            .warnings
            .iter()
            .any(|w| w.contains("basename differs"))
    );
    let prepared = export::prepare(&m, &Format::SnapGene, &review.context_revision, false).unwrap();
    let mut bad = review.clone();
    bad.context_revision = "0".repeat(64);
    assert!(bad.validate_binding(&review.context_revision).is_err());
    let mut bad = review.clone();
    bad.profile = "foreign".into();
    assert!(bad.validate_binding(&review.context_revision).is_err());
    let mut bad = review.clone();
    bad.filename = "../outside.dna".into();
    assert!(bad.validate_binding(&review.context_revision).is_err());
    let mut bad = prepared.clone();
    bad.bytes = vec![0; 128];
    assert!(bad.validate_binding(&review.context_revision).is_err());
    let mut bad = prepared;
    bad.bytes.clear();
    assert!(bad.validate_binding(&review.context_revision).is_err());
    let mut raw = serde_json::to_value(review).unwrap();
    raw["unexpected"] = "no authority".into();
    assert!(serde_json::from_value::<export::Assessment>(raw).is_err());
}
#[test]
fn worker_review_and_prepare_append_only_tags_use_the_same_context_and_old_header_refusal() {
    let m = molecule();
    let format = Format::Xdna;
    let context = export::revision(&m, &format).unwrap();
    let request = protocol::Request {
        id: 174,
        revision: context.clone(),
        job: protocol::AnalysisJob::ReviewSequenceExport {
            molecule: m.clone(),
            format: format.clone(),
        },
    };
    assert_eq!(
        postcard::take_from_bytes::<u32>(&postcard::to_allocvec(&request.job).unwrap())
            .unwrap()
            .0,
        173
    );
    let bytes = protocol::encode_request(&request).unwrap();
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&89u16.to_le_bytes());
    assert!(protocol::process_request(&old).is_err());
    let response = protocol::decode_response(&protocol::process_request(&bytes).unwrap()).unwrap();
    assert_eq!(
        postcard::take_from_bytes::<u32>(&postcard::to_allocvec(&response.result).unwrap())
            .unwrap()
            .0,
        163
    );
    let protocol::AnalysisResult::SequenceExportAssessment(review) = response.result else {
        panic!("export assessment")
    };
    review.validate_binding(&context).unwrap();
    let request = protocol::Request {
        id: 175,
        revision: context.clone(),
        job: protocol::AnalysisJob::PrepareSequenceExport {
            molecule: m,
            format,
            expected_context_revision: context.clone(),
            allow_omissions: true,
        },
    };
    assert_eq!(
        postcard::take_from_bytes::<u32>(&postcard::to_allocvec(&request.job).unwrap())
            .unwrap()
            .0,
        174
    );
    let response = protocol::process_request(&protocol::encode_request(&request).unwrap()).unwrap();
    let mut old = response.clone();
    old[4..6].copy_from_slice(&89u16.to_le_bytes());
    assert!(protocol::decode_response(&old).is_err());
    let response = protocol::decode_response(&response).unwrap();
    assert_eq!(
        postcard::take_from_bytes::<u32>(&postcard::to_allocvec(&response.result).unwrap())
            .unwrap()
            .0,
        164
    );
    let protocol::AnalysisResult::SequenceExportPrepared(prepared) = response.result else {
        panic!("prepared export")
    };
    prepared.validate_binding(&context).unwrap();
    same(&review, &prepared.assessment);
}
