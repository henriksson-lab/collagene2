use collagene_core::sequence::annotations::library::*;
use collagene_core::sequence::annotations::merge;
use collagene_model::{Alphabet, Feature, Location, Molecule, Sequence, Topology};

fn merge_request(
    target: &Molecule,
    source: merge::Source,
    selections: Vec<merge::Selection>,
) -> merge::Request {
    merge::Request {
        schema_version: 1,
        target_revision: target.revision(),
        source,
        selections,
        limits: Default::default(),
    }
}
fn source(name: &str, bases: &str, location: Location) -> Molecule {
    let mut molecule = Molecule::dna(name, bases, Topology::Linear).unwrap();
    molecule.features.push(Feature {
        kind: "misc_feature".into(),
        location,
        qualifiers: vec![("label".into(), Some("reference part".into()))],
        needs_review: false,
    });
    molecule
}
fn library(template: &str) -> FeatureLibrary {
    let mut library = FeatureLibrary::new("Independent parts", "2026.1").unwrap();
    library
        .add_feature(
            &source("source", template, Location::range(0, template.len())),
            0,
        )
        .unwrap();
    library
}
fn options() -> ScanOptions {
    ScanOptions {
        max_comparisons: 100_000,
        max_hits: 100,
        ..ScanOptions::default()
    }
}

#[test]
fn exact_joined_reverse_and_origin_sources_extract_in_biological_order_and_add_idempotently() {
    let mut source = source(
        "source",
        "AAAACCCCGGGGTTTT",
        Location::Complement(Box::new(Location::Join(vec![
            Location::range(12, 16),
            Location::range(0, 4),
        ]))),
    );
    source.topology = Topology::Circular;
    source.features[0].needs_review = true;
    let mut lib = FeatureLibrary::new("My library", "v1").unwrap();
    let id = lib.add_feature(&source, 0).unwrap();
    let revision = lib.revision();
    assert_eq!(lib.entries[0].sequence.as_str(), "TTTTAAAA");
    assert_eq!(lib.entries[0].source.location, source.features[0].location);
    assert_eq!(lib.entries[0].source.document_revision, source.revision());
    assert!(lib.entries[0].source.source_needs_review);
    assert_eq!(lib.add_feature(&source, 0).unwrap(), id);
    assert_eq!(lib.entries.len(), 1);
    assert_eq!(lib.revision(), revision);
    let bytes = save_library(&lib).unwrap();
    let restored = load_library(&bytes).unwrap();
    assert_eq!(save_library(&restored).unwrap(), bytes);
    assert_eq!(restored.revision(), revision);
}

#[test]
fn unsupported_source_locations_alphabets_and_bad_indices_are_atomic() {
    let mut lib = library("ACG");
    let before = lib.clone();
    let mut invalid = source(
        "bad",
        "ACGT",
        Location::Order(vec![Location::range(0, 2), Location::range(2, 4)]),
    );
    assert!(lib.add_feature(&invalid, 0).is_err());
    assert_eq!(lib, before);
    invalid.features[0].location = Location::Range {
        start: 0,
        end: 3,
        before: true,
        after: false,
    };
    assert!(lib.add_feature(&invalid, 0).is_err());
    assert_eq!(lib, before);
    invalid.features[0].location = Location::Unresolved("other:1..3".into());
    assert!(lib.add_feature(&invalid, 0).is_err());
    assert_eq!(lib, before);
    invalid.sequence = Sequence::new(Alphabet::Rna, "ACGU").unwrap();
    assert!(lib.add_feature(&invalid, 0).is_err());
    assert!(lib.add_feature(&before_source(), 1).is_err());
    assert_eq!(lib, before);
    fn before_source() -> Molecule {
        source("source", "ACG", Location::range(0, 3))
    }
}

#[test]
fn strict_json_rejects_tamper_duplicates_future_fields_and_nested_ignored_model_fields() {
    let lib = library("ACG");
    let raw = serde_json::to_value(&lib).unwrap();
    for mode in 0..8 {
        let mut changed = raw.clone();
        match mode {
            0 => changed["entries"][0]["sequence"]["bases"] = serde_json::json!("ACT"),
            1 => {
                changed["entries"][0]["source"]["document_revision"] =
                    serde_json::json!("0".repeat(64))
            }
            2 => changed["schema_version"] = serde_json::json!(2),
            3 => changed["future_inference"] = serde_json::json!(true),
            4 => changed["entries"][0]["sequence"]["hidden_chemistry"] = serde_json::json!("ssDNA"),
            5 => {
                changed["entries"][0]["source"]["location"]["Range"]["unrecorded_phase"] =
                    serde_json::json!(1)
            }
            6 => changed["entries"][0]["sequence"]["bases"] = serde_json::json!("acg"),
            _ => changed["entries"]
                .as_array_mut()
                .unwrap()
                .push(raw["entries"][0].clone()),
        }
        assert!(
            load_library(&serde_json::to_vec(&changed).unwrap()).is_err(),
            "mode {mode}"
        );
    }
    let bytes = String::from_utf8(save_library(&lib).unwrap()).unwrap();
    let duplicate = bytes.replacen(
        "\"schema_version\":1",
        "\"schema_version\":1,\"schema_version\":1",
        1,
    );
    assert!(load_library(duplicate.as_bytes()).is_err());
    let mut whitespace = b" \n".to_vec();
    whitespace.extend(save_library(&lib).unwrap());
    assert_eq!(load_library(&whitespace).unwrap(), lib);
    assert!(load_library(&vec![b' '; MAX_LIBRARY_BYTES + 1]).is_err());
}

#[test]
fn all_circular_rotations_and_reverse_hits_have_independent_expected_geometry() {
    let lib = library("CGTA");
    let reference = "CGTAAC";
    for rotation in 0..reference.len() {
        let rotated = format!("{}{}", &reference[rotation..], &reference[..rotation]);
        for reverse_reference in [false, true] {
            let sequence = if reverse_reference {
                independent_reverse(&rotated)
            } else {
                rotated.clone()
            };
            let target = Molecule::dna("circular", &sequence, Topology::Circular).unwrap();
            let report = scan(&target, &lib, &options()).unwrap();
            let expected = oracle(&sequence, "CGTA", true, true, 0);
            assert_eq!(
                tuples(&report),
                expected,
                "rotation {rotation}, reverse {reverse_reference}"
            );
            assert!(report.complete);
            assert_eq!(report.hits.len(), 1);
            for hit in report.hits {
                assert_eq!(hit.location.extract(&target.sequence).unwrap(), "CGTA");
                assert_eq!(hit.feature.location, hit.location);
                assert!(hit.feature.needs_review);
                assert_eq!(hit.coverage, 1.);
                assert_eq!(hit.identity_lower, 1.);
                assert_eq!(hit.identity_upper, 1.);
            }
        }
    }
}

// Independent oracle uses literal allele sets and forward slicing, not the scanner's
// bit masks, location builder, matching utility or reverse-complement implementation.
fn alleles(base: u8) -> &'static [u8] {
    match base {
        b'A' => b"A",
        b'C' => b"C",
        b'G' => b"G",
        b'T' => b"T",
        b'R' => b"AG",
        b'Y' => b"CT",
        b'S' => b"CG",
        b'W' => b"AT",
        b'K' => b"GT",
        b'M' => b"AC",
        b'B' => b"CGT",
        b'D' => b"AGT",
        b'H' => b"ACT",
        b'V' => b"ACG",
        b'N' => b"ACGT",
        _ => panic!("oracle alphabet"),
    }
}
fn independent_reverse(text: &str) -> String {
    text.bytes()
        .rev()
        .map(|base| match base {
            b'A' => 'T',
            b'T' => 'A',
            b'C' => 'G',
            b'G' => 'C',
            b'R' => 'Y',
            b'Y' => 'R',
            b'S' => 'S',
            b'W' => 'W',
            b'K' => 'M',
            b'M' => 'K',
            b'B' => 'V',
            b'D' => 'H',
            b'H' => 'D',
            b'V' => 'B',
            b'N' => 'N',
            _ => panic!("oracle complement"),
        })
        .collect()
}
fn expanded(base: u8) -> &'static [u8] {
    if base == b'Y' { b"CT" } else { alleles(base) }
}
fn oracle(
    reference: &str,
    query: &str,
    circular: bool,
    both: bool,
    max_mismatch: usize,
) -> Vec<(usize, bool, usize, usize, usize)> {
    let mut out = vec![];
    if query.len() > reference.len() {
        return out;
    }
    let window = if circular {
        format!("{reference}{reference}")
    } else {
        reference.into()
    };
    let starts = if circular {
        reference.len()
    } else {
        reference.len() - query.len() + 1
    };
    let reverse = independent_reverse(query);
    for start in 0..starts {
        for (is_reverse, template) in [(false, query), (true, reverse.as_str())]
            .into_iter()
            .take(if both { 2 } else { 1 })
        {
            let (mut concrete, mut ambiguous, mut mismatch) = (0, 0, 0);
            for (a, b) in template
                .bytes()
                .zip(window[start..start + query.len()].bytes())
            {
                let (a, b) = (expanded(a), expanded(b));
                if !a.iter().any(|base| b.contains(base)) {
                    mismatch += 1;
                } else if a.len() == 1 && b.len() == 1 {
                    concrete += 1;
                } else {
                    ambiguous += 1;
                }
            }
            if mismatch <= max_mismatch {
                out.push((start, is_reverse, concrete, ambiguous, mismatch));
            }
        }
    }
    out
}
fn tuples(report: &ScanReport) -> Vec<(usize, bool, usize, usize, usize)> {
    report
        .hits
        .iter()
        .map(|hit| {
            (
                hit.start,
                hit.reverse,
                hit.concrete_matches,
                hit.compatible_ambiguous,
                hit.mismatches,
            )
        })
        .collect()
}
fn words(alphabet: &[u8], length: usize) -> Vec<String> {
    let mut words = vec![String::new()];
    for _ in 0..length {
        words = words
            .into_iter()
            .flat_map(|word| {
                alphabet
                    .iter()
                    .map(move |base| format!("{word}{}", *base as char))
            })
            .collect();
    }
    words
}
#[test]
fn exhaustive_4096_small_ambiguity_orientation_topology_and_mismatch_oracle_cases() {
    let libraries = words(b"ACRN", 2)
        .into_iter()
        .map(|query| {
            let lib = library(&query);
            (query, lib)
        })
        .collect::<Vec<_>>();
    let mut cases = 0;
    for reference in words(b"ACRN", 3) {
        for (query, lib) in &libraries {
            for circular in [false, true] {
                for max_mismatch in [0, 1] {
                    let target = Molecule::dna(
                        "oracle",
                        &reference,
                        if circular {
                            Topology::Circular
                        } else {
                            Topology::Linear
                        },
                    )
                    .unwrap();
                    let mut opts = options();
                    opts.max_mismatches = max_mismatch;
                    let report = scan(&target, lib, &opts).unwrap();
                    assert_eq!(
                        tuples(&report),
                        oracle(&reference, query, circular, true, max_mismatch),
                        "{reference} / {query} / circle {circular} / mm {max_mismatch}"
                    );
                    assert!(report.complete);
                    for hit in &report.hits {
                        assert_eq!(hit.length, 2);
                        assert_eq!(
                            hit.concrete_matches + hit.compatible_ambiguous + hit.mismatches,
                            2
                        );
                        assert_eq!(hit.identity_lower, hit.concrete_matches as f64 / 2.);
                        assert_eq!(
                            hit.identity_upper,
                            (hit.concrete_matches + hit.compatible_ambiguous) as f64 / 2.
                        );
                        assert!(hit.feature.needs_review);
                    }
                    cases += 1;
                }
            }
        }
    }
    assert_eq!(cases, 4096);
}

#[test]
fn palindromes_and_repeated_locations_remain_separate_review_required_proposals() {
    let lib = library("AT");
    let target = Molecule::dna("repeat", "ATAT", Topology::Linear).unwrap();
    let report = scan(&target, &lib, &options()).unwrap();
    assert_eq!(
        tuples(&report),
        vec![
            (0, false, 2, 0, 0),
            (0, true, 2, 0, 0),
            (2, false, 2, 0, 0),
            (2, true, 2, 0, 0)
        ]
    );
    let ids = report
        .hits
        .iter()
        .map(|hit| hit.id.clone())
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(ids.len(), 4);
    assert!(report.hits.iter().all(|hit| {
        hit.warnings
            .iter()
            .any(|w| w.contains("Multiple locations"))
    }));
    assert!(report.hits.iter().all(|hit| {
        hit.warnings
            .iter()
            .any(|w| w.contains("own reverse complement"))
    }));
    let ambiguous = scan(
        &Molecule::dna("uncertain", "AN", Topology::Linear).unwrap(),
        &library("AR"),
        &ScanOptions {
            both_strands: false,
            ..options()
        },
    )
    .unwrap();
    assert_eq!(ambiguous.hits[0].identity_lower, 0.5);
    assert_eq!(ambiguous.hits[0].identity_upper, 1.);
    assert_eq!(ambiguous.hits[0].compatible_ambiguous, 1);
}

#[test]
fn search_limits_never_promote_partial_or_unsearched_placements() {
    let target = Molecule::dna("target", "ACG", Topology::Linear).unwrap();
    let lib = library("ACG");
    let partial = scan(
        &target,
        &lib,
        &ScanOptions {
            both_strands: false,
            max_comparisons: 2,
            ..options()
        },
    )
    .unwrap();
    assert!(!partial.complete);
    assert_eq!(partial.comparisons, 2);
    assert!(partial.hits.is_empty());
    let exact = scan(
        &target,
        &lib,
        &ScanOptions {
            both_strands: false,
            max_comparisons: 3,
            max_hits: 1,
            ..options()
        },
    )
    .unwrap();
    assert!(exact.complete);
    assert_eq!(exact.hits.len(), 1);
    assert_eq!(exact.comparisons, 3);
    let repeated = scan(
        &Molecule::dna("target", "AAA", Topology::Linear).unwrap(),
        &library("A"),
        &ScanOptions {
            both_strands: false,
            max_hits: 1,
            ..options()
        },
    )
    .unwrap();
    assert!(!repeated.complete);
    assert_eq!(repeated.hits.len(), 1);
    assert_eq!(repeated.comparisons, 2);
    let overlong = scan(&target, &library("ACGT"), &options()).unwrap();
    assert!(overlong.complete);
    assert!(overlong.hits.is_empty());
    assert_eq!(overlong.comparisons, 0);
    for opts in [
        ScanOptions {
            max_comparisons: 0,
            ..options()
        },
        ScanOptions {
            max_comparisons: MAX_SCAN_COMPARISONS + 1,
            ..options()
        },
        ScanOptions {
            max_hits: 0,
            ..options()
        },
        ScanOptions {
            max_hits: MAX_SCAN_HITS + 1,
            ..options()
        },
        ScanOptions {
            max_mismatches: 1001,
            ..options()
        },
    ] {
        assert!(scan(&target, &lib, &opts).is_err());
    }
}

#[test]
fn transferred_annotations_keep_source_evidence_but_strip_phase_translation_ids_and_unknown_coordinates()
 {
    let mut original = source("source", "ACG", Location::range(0, 3));
    original.features[0].kind = "CDS".into();
    original.features[0].qualifiers.extend([
        ("translation".into(), Some("SOURCE_PROTEIN".into())),
        ("codon_start".into(), Some("2".into())),
        ("transl_table".into(), Some("11".into())),
        ("protein_id".into(), Some("SOURCE_ID".into())),
        ("transl_except".into(), Some("(pos:1..3,aa:Sec)".into())),
        ("future_coordinate_evidence".into(), Some("1..3".into())),
        (
            "note".into(),
            Some("User-created descriptive evidence".into()),
        ),
        ("ApEinfo_fwdcolor".into(), Some("#aabbcc".into())),
    ]);
    let original_before = original.clone();
    let mut lib = FeatureLibrary::new("User library", "1").unwrap();
    lib.add_feature(&original, 0).unwrap();
    assert_eq!(lib.entries[0].qualifiers, original.features[0].qualifiers);
    let target = Molecule::dna("target", "TTACGAA", Topology::Linear).unwrap();
    let report = scan(&target, &lib, &options()).unwrap();
    assert_eq!(report.hits.len(), 1);
    let hit = &report.hits[0];
    assert_eq!(hit.location, Location::range(2, 5));
    let edited = merge::apply(
        &target,
        &merge_request(
            &target,
            merge::Source::Library {
                library: (&lib).clone(),
                report: (&report).clone(),
            },
            (&[hit.id.clone()])
                .iter()
                .map(|id| merge::Selection {
                    proposal_id: id.clone(),
                    decision: merge::Decision::Append,
                })
                .collect(),
        ),
    )
    .map(|applied| applied.molecule)
    .unwrap();
    assert_eq!(edited.sequence, target.sequence);
    assert_eq!(edited.features.len(), 1);
    let feature = &edited.features[0];
    assert!(feature.needs_review);
    assert!(feature.qualifiers.iter().any(|(key, _)| key == "note"));
    assert!(
        feature
            .qualifiers
            .iter()
            .any(|(key, _)| key == "ApEinfo_fwdcolor")
    );
    for key in [
        "translation",
        "codon_start",
        "transl_table",
        "protein_id",
        "transl_except",
        "future_coordinate_evidence",
    ] {
        assert!(
            !feature.qualifiers.iter().any(|(actual, _)| actual == key),
            "{key}"
        );
    }
    let ledger =
        collagene_core::sequence::annotations::merge::Ledger::from_feature(feature).unwrap();
    assert_eq!(ledger.contributions.len(), 1);
    let collagene_core::sequence::annotations::merge::SourcePin::Library {
        library_revision,
        source_revision,
        source_feature_index,
        proposal_id,
        ..
    } = &ledger.contributions[0].source
    else {
        panic!("not library attribution")
    };
    assert_eq!(library_revision, &lib.revision());
    assert_eq!(source_revision, &original.revision());
    assert_eq!(*source_feature_index, 0);
    assert_eq!(proposal_id, &hit.id);
    assert_eq!(original, original_before);
    assert!(target.features.is_empty());
    let restored = collagene_core::protocol::load_document(
        &collagene_core::protocol::save_document(&edited).unwrap(),
    )
    .unwrap();
    assert_eq!(restored, edited);
}

#[test]
fn stale_or_forged_reports_cannot_authorize_annotations_and_rejections_are_atomic() {
    let lib = library("ACG");
    let target = Molecule::dna("target", "TTACGAA", Topology::Linear).unwrap();
    let report = scan(&target, &lib, &options()).unwrap();
    let id = report.hits[0].id.clone();
    let mut changed = target.clone();
    changed.name = "renamed".into();
    assert!(
        merge::apply(
            &changed,
            &merge_request(
                &changed,
                merge::Source::Library {
                    library: (&lib).clone(),
                    report: (&report).clone()
                },
                (&[id.clone()])
                    .iter()
                    .map(|id| merge::Selection {
                        proposal_id: id.clone(),
                        decision: merge::Decision::Append
                    })
                    .collect()
            )
        )
        .map(|applied| applied.molecule)
        .is_err()
    );
    let mut changed_lib = lib.clone();
    changed_lib.version = "2".into();
    assert!(
        merge::apply(
            &target,
            &merge_request(
                &target,
                merge::Source::Library {
                    library: (&changed_lib).clone(),
                    report: (&report).clone()
                },
                (&[id.clone()])
                    .iter()
                    .map(|id| merge::Selection {
                        proposal_id: id.clone(),
                        decision: merge::Decision::Append
                    })
                    .collect()
            )
        )
        .map(|applied| applied.molecule)
        .is_err()
    );
    for mode in 0..7 {
        let mut forged = report.clone();
        match mode {
            0 => forged.hits[0].feature.kind = "gene".into(),
            1 => forged.hits[0].feature.location = Location::range(0, 3),
            2 => forged.hits[0].identity_lower = 0.5,
            3 => forged.complete = false,
            4 => forged.hits[0].feature.needs_review = false,
            5 => forged.hits[0].id = "fake".into(),
            _ => forged.options.max_mismatches = 1,
        }
        assert!(
            merge::apply(
                &target,
                &merge_request(
                    &target,
                    merge::Source::Library {
                        library: (&lib).clone(),
                        report: (&forged).clone()
                    },
                    (&[id.clone()])
                        .iter()
                        .map(|id| merge::Selection {
                            proposal_id: id.clone(),
                            decision: merge::Decision::Append
                        })
                        .collect()
                )
            )
            .map(|applied| applied.molecule)
            .is_err(),
            "mode {mode}"
        );
        assert!(target.features.is_empty());
    }
    assert!(
        merge::apply(
            &target,
            &merge_request(
                &target,
                merge::Source::Library {
                    library: (&lib).clone(),
                    report: (&report).clone()
                },
                Vec::<merge::Selection>::new()
            )
        )
        .map(|applied| applied.molecule)
        .is_err()
    );
    assert!(
        merge::apply(
            &target,
            &merge_request(
                &target,
                merge::Source::Library {
                    library: (&lib).clone(),
                    report: (&report).clone()
                },
                (&[id.clone(), id.clone()])
                    .iter()
                    .map(|id| merge::Selection {
                        proposal_id: id.clone(),
                        decision: merge::Decision::Append
                    })
                    .collect()
            )
        )
        .map(|applied| applied.molecule)
        .is_err()
    );
    assert!(
        merge::apply(
            &target,
            &merge_request(
                &target,
                merge::Source::Library {
                    library: (&lib).clone(),
                    report: (&report).clone()
                },
                vec![merge::Selection {
                    proposal_id: "missing".into(),
                    decision: merge::Decision::Append,
                }]
            )
        )
        .map(|applied| applied.molecule)
        .is_err()
    );
    assert!(target.features.is_empty());
    // JSON itself does not make a legitimate report authoritative: reconstruction is still required.
    let reopened: ScanReport =
        serde_json::from_slice(&serde_json::to_vec(&report).unwrap()).unwrap();
    assert_eq!(
        merge::apply(
            &target,
            &merge_request(
                &target,
                merge::Source::Library {
                    library: (&lib).clone(),
                    report: (&reopened).clone()
                },
                (&[id])
                    .iter()
                    .map(|id| merge::Selection {
                        proposal_id: id.clone(),
                        decision: merge::Decision::Append
                    })
                    .collect()
            )
        )
        .map(|applied| applied.molecule)
        .unwrap()
        .features
        .len(),
        1
    );
}

#[test]
fn existing_and_cross_entry_duplicate_annotations_are_explicit_and_never_overwritten() {
    let mut lib = library("ACG");
    let mut another = source("other source", "ACG", Location::range(0, 3));
    another.features[0]
        .qualifiers
        .push(("note".into(), Some("other evidence".into())));
    lib.add_feature(&another, 0).unwrap();
    let target = Molecule::dna("target", "TTACGAA", Topology::Linear).unwrap();
    let report = scan(&target, &lib, &options()).unwrap();
    assert_eq!(report.hits.len(), 2);
    let ids = report
        .hits
        .iter()
        .map(|hit| hit.id.clone())
        .collect::<Vec<_>>();
    let merged = merge::apply(
        &target,
        &merge_request(
            &target,
            merge::Source::Library {
                library: lib.clone(),
                report: report.clone(),
            },
            ids.iter()
                .map(|id| merge::Selection {
                    proposal_id: id.clone(),
                    decision: merge::Decision::Append,
                })
                .collect(),
        ),
    )
    .unwrap()
    .molecule;
    assert_eq!(merged.features.len(), 2);
    assert!(
        merged
            .features
            .iter()
            .all(|feature| feature.location == Location::range(2, 5) && feature.needs_review)
    );
    assert_eq!(
        merged
            .features
            .iter()
            .filter(|feature| feature
                .qualifiers
                .iter()
                .any(|(key, value)| key == "note" && value.as_deref() == Some("other evidence")))
            .count(),
        1
    );
    let contributions = merged
        .features
        .iter()
        .map(|feature| merge::Ledger::from_feature(feature).unwrap())
        .collect::<Vec<_>>();
    assert!(
        contributions
            .iter()
            .all(|ledger| ledger.contributions.len() == 1)
    );
    assert_ne!(
        contributions[0].contributions[0].source,
        contributions[1].contributions[0].source
    );
    assert!(target.features.is_empty());
    let mut existing = target.clone();
    existing.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::Join(vec![Location::range(2, 3), Location::range(3, 5)]),
        qualifiers: vec![
            ("label".into(), Some("reference part".into())),
            ("note".into(), Some("Do not replace".into())),
        ],
        needs_review: false,
    });
    let before = existing.clone();
    let matched = scan(&existing, &lib, &options()).unwrap();
    assert!(matched.hits.iter().all(|hit| hit.same_label_geometry_hint));
    let extended = merge::apply(
        &existing,
        &merge_request(
            &existing,
            merge::Source::Library {
                library: lib.clone(),
                report: matched.clone(),
            },
            vec![merge::Selection {
                proposal_id: matched.hits[0].id.clone(),
                decision: merge::Decision::Append,
            }],
        ),
    )
    .unwrap()
    .molecule;
    assert_eq!(extended.features.len(), 2);
    assert_eq!(extended.features[0], before.features[0]);
    assert!(extended.features[1].needs_review);
    assert_eq!(
        merge::Ledger::from_feature(&extended.features[1])
            .unwrap()
            .contributions
            .len(),
        1
    );
    assert_eq!(existing, before);
}

#[test]
fn empty_libraries_and_invalid_library_fields_have_explicit_bounded_behavior() {
    let lib = FeatureLibrary::new("empty", "1").unwrap();
    let target = Molecule::dna("target", "ACG", Topology::Linear).unwrap();
    let report = scan(&target, &lib, &options()).unwrap();
    assert!(report.complete);
    assert!(report.hits.is_empty());
    assert!(report.warnings.iter().any(|w| w.contains("no entries")));
    assert!(FeatureLibrary::new("", "1").is_err());
    assert!(FeatureLibrary::new("library", "").is_err());
    let mut invalid = library("ACG");
    invalid.entries[0].revision = "0".repeat(64);
    assert!(invalid.validate().is_err());
    let mut invalid = library("ACG");
    invalid.entries[0].qualifiers.push((
        "note".into(),
        Some("X".repeat(MAX_QUALIFIER_VALUE_BYTES + 1)),
    ));
    assert!(invalid.validate().is_err());
    let protein = Molecule {
        name: "protein".into(),
        sequence: Sequence::new(Alphabet::Protein, "MKW").unwrap(),
        topology: Topology::Linear,
        features: vec![],
        metadata: Default::default(),
        original_genbank: None,
    };
    assert!(scan(&protein, &lib, &options()).is_err());
}

#[test]
fn long_cds_translation_remains_exact_library_evidence_without_being_copied_as_current_translation()
{
    let mut original = source("long evidence", "ACG", Location::range(0, 3));
    let translation = "M".repeat(16_384);
    original.features[0].kind = "CDS".into();
    original.features[0]
        .qualifiers
        .push(("translation".into(), Some(translation.clone())));
    let mut lib = FeatureLibrary::new("Large source qualifiers", "1").unwrap();
    lib.add_feature(&original, 0).unwrap();
    let reopened = load_library(&save_library(&lib).unwrap()).unwrap();
    assert_eq!(
        reopened.entries[0]
            .qualifiers
            .iter()
            .find(|(key, _)| key == "translation")
            .unwrap()
            .1
            .as_ref(),
        Some(&translation)
    );
    let target = Molecule::dna("target", "ACG", Topology::Linear).unwrap();
    let report = scan(&target, &reopened, &options()).unwrap();
    assert!(
        !report.hits[0]
            .feature
            .qualifiers
            .iter()
            .any(|(key, _)| key == "translation")
    );
    assert!(
        report.hits[0]
            .warnings
            .iter()
            .any(|warning| warning.contains("translation"))
    );
}

#[test]
fn ambiguous_identity_report_json_preserves_fraction_bits_and_remains_reconstructable() {
    let lib = library("ACG");
    let target = Molecule::dna("ambiguous", "ANG", Topology::Linear).unwrap();
    let report = scan(&target, &lib, &options()).unwrap();
    assert_eq!(report.hits.len(), 1);
    assert_eq!(report.hits[0].identity_lower, 2. / 3.);
    assert_eq!(report.hits[0].identity_upper, 1.);
    let bytes = serde_json::to_vec(&report).unwrap();
    let restored: ScanReport = serde_json::from_slice(&bytes).unwrap();
    assert_eq!(
        restored.hits[0].identity_lower.to_bits(),
        report.hits[0].identity_lower.to_bits()
    );
    assert_eq!(serde_json::to_vec(&restored).unwrap(), bytes);
    let result = merge::apply(
        &target,
        &merge_request(
            &target,
            merge::Source::Library {
                library: (&lib).clone(),
                report: (&restored).clone(),
            },
            (&[restored.hits[0].id.clone()])
                .iter()
                .map(|id| merge::Selection {
                    proposal_id: id.clone(),
                    decision: merge::Decision::Append,
                })
                .collect(),
        ),
    )
    .map(|applied| applied.molecule)
    .unwrap();
    assert!(result.features[0].needs_review);
    let ledger = merge::Ledger::from_feature(&result.features[0]).unwrap();
    assert_eq!(ledger.contributions.len(), 1);
    let merge::SourcePin::Library {
        concrete_matches,
        compatible_ambiguous,
        mismatches,
        template_length,
        ..
    } = &ledger.contributions[0].source
    else {
        panic!("not library evidence")
    };
    assert_eq!(
        (
            *concrete_matches,
            *compatible_ambiguous,
            *mismatches,
            *template_length
        ),
        (2, 1, 0, 3)
    );
    assert_eq!(
        *concrete_matches as f64 / *template_length as f64,
        report.hits[0].identity_lower
    );
}

#[test]
fn large_retained_descriptive_evidence_reaches_report_limit_without_oversized_output() {
    let mut original = source("large notes", "A", Location::range(0, 1));
    for _ in 0..4 {
        original.features[0]
            .qualifiers
            .push(("note".into(), Some("X".repeat(512 * 1024))));
    }
    let mut lib = FeatureLibrary::new("Large evidence", "1").unwrap();
    lib.add_feature(&original, 0).unwrap();
    let target = Molecule::dna("repeat", "AAAAAAAAAAAAAAAAAAAA", Topology::Linear).unwrap();
    let report = scan(
        &target,
        &lib,
        &ScanOptions {
            both_strands: false,
            ..options()
        },
    )
    .unwrap();
    assert!(!report.complete);
    assert!(!report.hits.is_empty());
    assert!(report.hits.len() < 20);
    assert!(serde_json::to_vec(&report).unwrap().len() <= 32 * 1024 * 1024);
    assert!(
        report
            .warnings
            .iter()
            .any(|warning| warning.contains("report limit"))
    );
    for hit in &report.hits {
        assert_eq!(hit.concrete_matches, 1);
        assert_eq!(hit.identity_lower, 1.);
        assert!(hit.feature.needs_review);
    }
}

#[test]
fn every_iupac_allele_pair_and_reverse_complement_matches_literal_set_oracle() {
    for query in b"ACGTRYSWKMBDHVN" {
        let query = (*query as char).to_string();
        let lib = library(&query);
        for reference in b"ACGTRYSWKMBDHVN" {
            let reference = (*reference as char).to_string();
            let target = Molecule::dna("allele pair", &reference, Topology::Linear).unwrap();
            for max_mismatch in [0, 1] {
                let opts = ScanOptions {
                    max_mismatches: max_mismatch,
                    ..options()
                };
                let report = scan(&target, &lib, &opts).unwrap();
                assert_eq!(
                    tuples(&report),
                    oracle(&reference, &query, false, true, max_mismatch),
                    "{query}/{reference}, mismatch budget {max_mismatch}"
                );
                assert!(report.complete);
            }
        }
    }
}
