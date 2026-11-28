//! Independent citation claims: bibliography survives without acquiring new coverage.
use collagene_model::genbank_edit::{
    CitationProjection, EditLog, KEY, SourceOnlyReason, VerifiedLog,
};
use collagene_model::{Molecule, Sequence, Topology};

fn source(description: &str) -> Molecule {
    let mut source = Molecule::dna("reference", "AAAACCCCGGGG", Topology::Linear).unwrap();
    source.original_genbank = Some(format!(
        "LOCUS       reference                 12 bp    DNA     linear   SYN 01-JAN-2026\nREFERENCE   {description}\n  AUTHORS   Example,A.\n  TITLE     Literal citation\n  JOURNAL   Unpublished\nORIGIN\n        1 aaaaccccgg gg\n//\n"
    ));
    source
}
fn verify(molecule: &Molecule) -> VerifiedLog {
    VerifiedLog::verify(molecule, &Sequence::dna("AAAACCCCGGGG").unwrap(), 0)
        .unwrap()
        .unwrap()
}
fn source_only(log: &VerifiedLog, description: &str, expected: SourceOnlyReason) {
    match log.project_references(&[description]).unwrap().remove(0) {
        CitationProjection::SourceOnly {
            number,
            original_description,
            reason,
        } => {
            assert_eq!(
                number,
                description
                    .split_whitespace()
                    .next()
                    .unwrap()
                    .parse::<u64>()
                    .unwrap()
            );
            assert_eq!(original_description, description);
            assert_eq!(reason, expected);
        }
        other => panic!("unexpected destination claim: {other:?}"),
    }
}

#[test]
fn legal_noncoordinate_scopes_never_become_destination_claims() {
    for (description, reason) in [
        ("1", SourceOnlyReason::Unspecified),
        ("1  (sites)", SourceOnlyReason::Sites),
        ("4294967296  (sites)", SourceOnlyReason::Sites),
        (
            "1  (bases 184467440737095516160 to 999999999999999999999; annotation prose)",
            SourceOnlyReason::UnsupportedScope,
        ),
        (
            "1  (annotation supplied by the submitter)",
            SourceOnlyReason::UnsupportedScope,
        ),
        (
            "1  (bases 3 to 8) trailing junk",
            SourceOnlyReason::UnsupportedScope,
        ),
        (
            "1  (bases 3 to 8; bases 10 to 12) trailing junk",
            SourceOnlyReason::UnsupportedScope,
        ),
        (
            "1  (bases 3 to 8; unexplained)",
            SourceOnlyReason::UnsupportedScope,
        ),
    ] {
        let original = source(description);
        let edited = original.replace(0, 0, "TT").unwrap();
        assert_eq!(edited.original_genbank, original.original_genbank);
        source_only(&verify(&edited), description, reason);
    }
}

#[test]
fn wholly_deleted_coverage_cannot_reappear_after_later_insertion() {
    // Newly inserted GG survives, but both original cited bases are deleted.
    let interrupted = source("1  (bases 3 to 4)")
        .replace(3, 3, "GG")
        .unwrap()
        .replace(2, 3, "")
        .unwrap()
        .replace(4, 5, "")
        .unwrap();
    assert_eq!(interrupted.sequence.as_str(), "AAGGCCCCGGGG");
    source_only(
        &verify(&interrupted),
        "1  (bases 3 to 4)",
        SourceOnlyReason::Deleted,
    );
    let replaced = source("1  (bases 3 to 4)").replace(2, 4, "TT").unwrap();
    source_only(
        &verify(&replaced),
        "1  (bases 3 to 4)",
        SourceOnlyReason::Deleted,
    );
    // Inserting and then removing a gap restores contiguous surviving bases.
    let restored = source("1  (bases 3 to 4)")
        .replace(3, 3, "GG")
        .unwrap()
        .replace(3, 5, "")
        .unwrap();
    assert_eq!(
        verify(&restored)
            .project_references(&["1  (bases 3 to 4)"])
            .unwrap(),
        vec![CitationProjection::Mapped {
            number: 1,
            description: "1  (bases 3 to 4)".into()
        }]
    );
    let deleted = source("1  (bases 3 to 4)").replace(2, 4, "").unwrap();
    source_only(
        &verify(&deleted),
        "1  (bases 3 to 4)",
        SourceOnlyReason::Deleted,
    );
    let inserted = deleted.replace(2, 2, "TT").unwrap();
    assert_eq!(inserted.sequence.as_str(), "AATTCCCCGGGG");
    source_only(
        &verify(&inserted),
        "1  (bases 3 to 4)",
        SourceOnlyReason::Deleted,
    );
}

#[test]
fn repeated_bases_keyword_has_literal_full_scope() {
    let edited = source("1  (bases 3 to 8; bases 10 to 12)")
        .replace(0, 0, "TT")
        .unwrap();
    assert_eq!(
        verify(&edited)
            .project_references(&["1  (bases 3 to 8; bases 10 to 12)"])
            .unwrap(),
        vec![CitationProjection::Mapped {
            number: 1,
            description: "1  (bases 5 to 10; 12 to 14)".into()
        }]
    );
}

#[test]
fn source_only_is_not_a_fallback_for_invalid_replay_authority() {
    let original = source("1  (sites)");
    let edited = original.replace(0, 0, "TT").unwrap();
    let baseline = &original.sequence;
    for mutation in 0..6 {
        let mut forged = edited.clone();
        let mut log = EditLog::read(&forged).unwrap().unwrap();
        match mutation {
            0 => {
                forged.original_genbank = Some(
                    original
                        .original_genbank
                        .as_ref()
                        .unwrap()
                        .replace("Literal", "Changed"),
                )
            }
            1 => forged.sequence = Sequence::dna("AAAAAACCCCGGGG").unwrap(),
            2 => log.version = 2,
            3 => log.baseline_sequence_blake3 = "0".repeat(64),
            4 => log.edits[0].replacement = "AA".into(),
            5 => log.source_record_index = 1,
            _ => unreachable!(),
        }
        forged
            .metadata
            .insert(KEY.into(), serde_json::to_string(&log).unwrap());
        assert!(
            VerifiedLog::verify(&forged, baseline, 0).is_err(),
            "mutation {mutation}"
        );
    }
    assert!(VerifiedLog::verify(&edited, &Sequence::dna("CCCCAAAAGGGG").unwrap(), 0).is_err());
    assert!(VerifiedLog::verify(&edited, baseline, 1).is_err());
    let mut aligned_forgery = edited.clone();
    let mut log = EditLog::read(&aligned_forgery).unwrap().unwrap();
    log.source_record_index = 1;
    aligned_forgery
        .metadata
        .insert(KEY.into(), serde_json::to_string(&log).unwrap());
    aligned_forgery
        .metadata
        .insert("source_record_index".into(), "1".into());
    assert!(VerifiedLog::verify(&aligned_forgery, baseline, 1).is_err());
    assert!(aligned_forgery.replace(0, 0, "A").is_err());
}

#[test]
fn citation_bounds_apply_to_source_only_and_mapped_descriptions() {
    let edited = source("1").replace(0, 0, "TT").unwrap();
    let log = verify(&edited);
    let bare: Vec<String> = (1..=1024).map(|number| number.to_string()).collect();
    assert_eq!(log.project_references(&bare).unwrap().len(), 1024);
    assert!(
        log.project_references((1..=1025).map(|number| number.to_string()))
            .is_err()
    );
    let max_description = format!("1 {}", "x".repeat(4094));
    assert!(log.project_references(&[max_description.as_str()]).is_ok());
    assert!(
        log.project_references(&[format!("{max_description}x").as_str()])
            .is_err()
    );
    let description = |n| format!("1  (bases {})", vec!["1 to 1"; n].join("; "));
    let spans = description(256);
    assert!(log.project_references(&[spans.as_str()]).is_ok());
    assert!(
        log.project_references(&[description(257).as_str()])
            .is_err()
    );
    let numbered_spans = |count| {
        (1..=count)
            .map(|number| format!("{number}  (bases {})", vec!["1 to 1"; 256].join("; ")))
            .collect::<Vec<_>>()
    };
    assert!(log.project_references(numbered_spans(16)).is_ok());
    assert!(log.project_references(numbered_spans(17)).is_err());
}

#[test]
fn citation_projection_work_is_bounded_independently_of_ledger_byte_size() {
    let mut edited = source("1");
    // Substitutions cannot coalesce into a typing insertion. 244 * 4096
    // projected spans is below one million; the next edit crosses that bound.
    for index in 0..244 {
        edited = edited
            .replace(10, 11, if index % 2 == 0 { "A" } else { "G" })
            .unwrap();
    }
    let descriptions: Vec<String> = (1..=16)
        .map(|number| format!("{number}  (bases {})", vec!["1 to 1"; 256].join("; ")))
        .collect();
    assert!(verify(&edited).project_references(&descriptions).is_ok());
    edited = edited.replace(10, 11, "A").unwrap();
    assert!(verify(&edited).project_references(&descriptions).is_err());
}

#[test]
fn numeric_scope_errors_and_invalid_reference_framing_are_rejected() {
    let edited = source("1").replace(0, 0, "TT").unwrap();
    let log = verify(&edited);
    for description in [
        "0  (sites)",
        "not-a-number  (sites)",
        "1\u{0001} (sites)",
        "1  (bases 0 to 3)",
        "1  (bases 8 to 3)",
        "1  (bases 1 to 13)",
    ] {
        assert!(
            log.project_references(&[description]).is_err(),
            "{description:?}"
        );
    }
}

#[test]
fn exact_source_coverage_representation_limit_retains_citation_as_source_only() {
    let baseline = Sequence::dna(&"A".repeat(600)).unwrap();
    let mut edited = Molecule::dna("fragmented", baseline.as_str(), Topology::Linear).unwrap();
    edited.original_genbank = Some(format!(
        "LOCUS       fragmented               600 bp    DNA     linear   SYN 01-JAN-2026\nREFERENCE   1  (bases 1 to 600)\nORIGIN\n        1 {}\n//\n",
        "a".repeat(600)
    ));
    // Distinct internal insertions produce 257 exact ancestral spans. Inserted
    // C bases must never be swallowed to fit a coordinate envelope.
    for index in 0..256 {
        edited = edited.replace(index * 2 + 1, index * 2 + 1, "C").unwrap();
    }
    let log = VerifiedLog::verify(&edited, &baseline, 0).unwrap().unwrap();
    source_only(
        &log,
        "1  (bases 1 to 600)",
        SourceOnlyReason::RepresentationLimit,
    );
}
