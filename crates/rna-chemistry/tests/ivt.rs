//! Independent literal coordinate and requested-constraint fixtures. These tests
//! do not equate sequence conformance with successful transcription or capping.
use collagene_model::{Alphabet, Molecule, Sequence, Topology};
use collagene_rna_chemistry::{Provenance, ivt::*};

fn evidence() -> Provenance {
    Provenance {
        statement: "Authored DNA geometry and intended transcript; not a measured material state"
            .into(),
        reference: None,
    }
}
fn source(dna: &str, topology: Topology) -> Molecule {
    Molecule::dna("Independent IVT template fixture", dna, topology).unwrap()
}
fn request(source: &Molecule, start: usize, length: usize) -> Request {
    Request {
        source_revision: source.revision(),
        profile: profiles().unwrap().remove(0).reference,
        transcript: Transcript {
            start,
            length,
            coding_strand: CodingStrand::Forward,
        },
        cap: CapChoice::CleanCapAgN7113,
        initiator: InitiatorRequirement::NoAdditionalRequirement,
        linearization: Some(Linearization {
            coding_cut_offset: length,
            template_cut_offset: length,
            provenance: evidence(),
        }),
        poly_a: None,
        provenance: evidence(),
    }
}
fn verdict(report: &Report, rule: Rule) -> Verdict {
    report
        .findings
        .iter()
        .find(|finding| finding.rule == rule)
        .unwrap()
        .verdict
}
fn fixture() -> (Molecule, Request) {
    // 2 upstream +17 promoter +13 transcript +2 downstream =34 bases.
    // Coding transcript AGGC AAAA G AAAA, poly(A) offsets4..13.
    let m = source("TTTAATACGACTCACTATAAGGCAAAAGAAAACC", Topology::Linear);
    let mut r = request(&m, 19, 13);
    r.poly_a = Some(PolyARegion {
        start_offset: 4,
        end_offset: 13,
    });
    (m, r)
}

#[test]
fn literal_ag_template_checks_and_encoded_region_keep_exact_source_coordinates() {
    let (m, r) = fixture();
    assert_eq!(m.sequence.len(), 34);
    let report = assess(&m, &r).unwrap();
    assert_eq!(report.promoter.dna, "TAATACGACTCACTATA");
    assert_eq!(
        report.promoter.original_positions,
        (2..19).collect::<Vec<_>>()
    );
    assert_eq!(report.initiator_context.dna, "AGGCAAAA");
    assert_eq!(
        report.transcript_spans,
        vec![MappedSpan {
            start: 19,
            end: 32,
            reverse: false
        }]
    );
    for rule in [
        Rule::ConsensusPromoter,
        Rule::CapInitiator,
        Rule::Linearization,
        Rule::RunoffBoundary,
        Rule::PolyATerminus,
    ] {
        assert_eq!(verdict(&report, rule), Verdict::Pass);
    }
    assert_eq!(
        verdict(&report, Rule::PolyAEncoding),
        Verdict::Informational
    );
    let encoding = report.poly_a.as_ref().unwrap();
    assert_eq!(
        (
            encoding.adenosines,
            encoding.other_concrete_bases,
            encoding.ambiguous_bases
        ),
        (8, 1, 0)
    );
    assert_eq!(encoding.longest_exact_a_run, 4);
    assert_eq!(encoding.runs_total, 3);
    assert!(encoding.segmentation_complete);
    assert_eq!(
        encoding
            .runs
            .iter()
            .map(|run| (run.start_offset, run.end_offset, run.dna_preview.as_str()))
            .collect::<Vec<_>>(),
        vec![(4, 8, "AAAA"), (8, 9, "G"), (9, 13, "AAAA")]
    );
    assert_eq!(
        encoding.runs[0].source_spans,
        vec![MappedSpan {
            start: 23,
            end: 27,
            reverse: false
        }]
    );
    assert_eq!(report.linearization.as_ref().unwrap().template_cut_bond, 32);
    report.validate_for(&m, &r).unwrap();
}

#[test]
fn independently_requested_ggg_is_a_conflict_with_ag_and_never_a_default() {
    let (m, mut r) = fixture();
    assert!(
        !assess(&m, &r)
            .unwrap()
            .findings
            .iter()
            .any(|finding| finding.rule == Rule::DeclaredInitiator)
    );
    r.initiator = InitiatorRequirement::Ggg;
    let ag = assess(&m, &r).unwrap();
    assert_eq!(verdict(&ag, Rule::CapInitiator), Verdict::Pass);
    assert_eq!(verdict(&ag, Rule::DeclaredInitiator), Verdict::Fail);
    assert_eq!(verdict(&ag, Rule::ConflictingRequirements), Verdict::Fail);
    assert_eq!(ag.initiator_context.dna, "AGGCAAAA");
    // Literal alternative source meets GGG but cannot meet +1A simultaneously.
    let g = source("TTTAATACGACTCACTATAGGGCAAAAGAAAACC", Topology::Linear);
    r.source_revision = g.revision();
    let report = assess(&g, &r).unwrap();
    assert_eq!(verdict(&report, Rule::CapInitiator), Verdict::Fail);
    assert_eq!(verdict(&report, Rule::DeclaredInitiator), Verdict::Pass);
    assert_eq!(
        verdict(&report, Rule::ConflictingRequirements),
        Verdict::Fail
    );
    r.cap = CapChoice::EnzymaticPostTranscriptional;
    let report = assess(&g, &r).unwrap();
    assert_eq!(verdict(&report, Rule::CapInitiator), Verdict::Informational);
    assert!(
        !report
            .findings
            .iter()
            .any(|finding| finding.rule == Rule::ConflictingRequirements)
    );
}

fn independent_reverse_complement(dna: &str) -> String {
    dna.bytes()
        .rev()
        .map(|base| match base {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            _ => panic!("concrete literal only"),
        })
        .collect()
}

#[test]
fn both_strands_and_every_circular_origin_preserve_literal_calls_spans_and_bonds() {
    let (forward, _) = fixture();
    let n = forward.sequence.len();
    for reverse in [false, true] {
        let original = if reverse {
            independent_reverse_complement(forward.sequence.as_str())
        } else {
            forward.sequence.as_str().to_owned()
        };
        let start = if reverse { 14 } else { 19 };
        for rotation in 0..n {
            let dna = format!("{}{}", &original[rotation..], &original[..rotation]);
            let m = source(&dna, Topology::Circular);
            let mut r = request(&m, (start + n - rotation) % n, 13);
            r.transcript.coding_strand = if reverse {
                CodingStrand::Reverse
            } else {
                CodingStrand::Forward
            };
            r.poly_a = Some(PolyARegion {
                start_offset: 4,
                end_offset: 13,
            });
            let report = assess(&m, &r).unwrap();
            assert_eq!(report.promoter.dna, "TAATACGACTCACTATA");
            assert_eq!(report.initiator_context.dna, "AGGCAAAA");
            let expected_positions = (0..13)
                .map(|offset| {
                    let before_rotation = if reverse { 14 - offset } else { 19 + offset };
                    (before_rotation + n - rotation) % n
                })
                .collect::<Vec<_>>();
            let mapped_positions = report
                .transcript_spans
                .iter()
                .flat_map(|span| {
                    let mut positions = (span.start..span.end).collect::<Vec<_>>();
                    if span.reverse {
                        positions.reverse();
                    }
                    positions
                })
                .collect::<Vec<_>>();
            assert_eq!(mapped_positions, expected_positions);
            let original_bond = if reverse { 2 } else { 32 };
            assert_eq!(
                report.linearization.unwrap().template_cut_bond,
                (original_bond + n - rotation) % n
            );
            assert_eq!(report.poly_a.unwrap().adenosines, 8);
        }
    }
}

#[test]
fn supplied_strand_cut_offsets_distinguish_five_prime_three_prime_and_runoff_mismatch() {
    let (m, mut r) = fixture();
    r.linearization.as_mut().unwrap().coding_cut_offset = 11;
    let report = assess(&m, &r).unwrap();
    let end = report.linearization.as_ref().unwrap();
    assert_eq!(
        (
            end.kind,
            end.overhang_bases,
            end.coding_cut_bond,
            end.template_cut_bond
        ),
        (EndKind::FivePrimeOverhang, 2, 30, 32)
    );
    assert_eq!(verdict(&report, Rule::Linearization), Verdict::Pass);
    r.linearization.as_mut().unwrap().coding_cut_offset = 15;
    let report = assess(&m, &r).unwrap();
    assert_eq!(
        report.linearization.as_ref().unwrap().kind,
        EndKind::ThreePrimeOverhang
    );
    assert_eq!(verdict(&report, Rule::Linearization), Verdict::Advisory);
    r.linearization.as_mut().unwrap().template_cut_offset = 12;
    assert_eq!(
        verdict(&assess(&m, &r).unwrap(), Rule::RunoffBoundary),
        Verdict::Fail
    );
    r.linearization = None;
    let report = assess(&m, &r).unwrap();
    assert!(report.linearization.is_none());
    assert_eq!(
        verdict(&report, Rule::Linearization),
        Verdict::Indeterminate
    );
}

#[test]
fn missing_and_ambiguous_bases_remain_indeterminate_without_invented_absence() {
    let m = source("AGGAAA", Topology::Linear);
    let mut r = request(&m, 0, 6);
    r.cap = CapChoice::NotSpecified;
    let report = assess(&m, &r).unwrap();
    assert_eq!(
        verdict(&report, Rule::ConsensusPromoter),
        Verdict::Indeterminate
    );
    assert_eq!(verdict(&report, Rule::CapInitiator), Verdict::Indeterminate);
    assert!(!report.promoter.complete);
    let m = source("TAATACGACTCACTATARNAAANAAA", Topology::Linear);
    r = request(&m, 17, 8);
    r.poly_a = Some(PolyARegion {
        start_offset: 2,
        end_offset: 8,
    });
    let report = assess(&m, &r).unwrap();
    assert_eq!(verdict(&report, Rule::CapInitiator), Verdict::Indeterminate);
    assert_eq!(
        verdict(&report, Rule::PolyAEncoding),
        Verdict::Indeterminate
    );
    assert_eq!(report.poly_a.unwrap().ambiguous_bases, 1);
}

#[test]
fn encoded_runs_are_bounded_but_full_counts_and_incompleteness_remain_explicit() {
    let m = source(
        &format!("TAATACGACTCACTATAAG{}", "AC".repeat(1100)),
        Topology::Linear,
    );
    let mut r = request(&m, 17, 2202);
    r.poly_a = Some(PolyARegion {
        start_offset: 2,
        end_offset: 2202,
    });
    let report = assess(&m, &r).unwrap();
    let tail = report.poly_a.as_ref().unwrap();
    assert_eq!(
        (
            tail.adenosines,
            tail.other_concrete_bases,
            tail.longest_exact_a_run
        ),
        (1100, 1100, 1)
    );
    assert_eq!(tail.runs_total, 2200);
    assert_eq!(tail.runs.len(), MAX_RETAINED_RUNS);
    assert!(!tail.segmentation_complete);
    report.validate_for(&m, &r).unwrap();
}

#[test]
fn nonterminal_encoded_region_is_not_silently_extended_to_the_transcript_end() {
    let (m, mut r) = fixture();
    r.poly_a.as_mut().unwrap().end_offset = 8;
    let report = assess(&m, &r).unwrap();
    assert_eq!(verdict(&report, Rule::PolyATerminus), Verdict::Fail);
    assert_eq!(report.poly_a.as_ref().unwrap().adenosines, 4);
    assert!(!report.poly_a.unwrap().terminal);
}

#[test]
fn immutable_source_profile_options_and_documentary_report_are_all_bound() {
    let (m, r) = fixture();
    let report = assess(&m, &r).unwrap();
    let reopened: Report = serde_json::from_slice(&serde_json::to_vec(&report).unwrap()).unwrap();
    reopened.validate_for(&m, &r).unwrap();
    let mut changed = r.clone();
    changed.source_revision = "0".repeat(64);
    assert!(assess(&m, &changed).is_err());
    changed = r.clone();
    changed.profile.identity = "0".repeat(64);
    assert!(assess(&m, &changed).is_err());
    changed = r.clone();
    changed.cap = CapChoice::Uncapped;
    assert!(report.validate_for(&m, &changed).is_err());
    let mut edited = report.clone();
    edited.profile.sources[0].url.push_str("/claimed");
    assert!(edited.validate_binding(&report.context_revision).is_err());
    edited = report.clone();
    edited.poly_a.as_mut().unwrap().adenosines += 1;
    assert!(edited.validate_for(&m, &r).is_err());
    let mut wire = serde_json::to_value(&r).unwrap();
    wire.as_object_mut().unwrap().remove("linearization");
    assert!(serde_json::from_value::<Request>(wire).is_err());
}

#[test]
fn resource_geometry_and_non_dna_inputs_are_refused_before_scientific_checks() {
    let (m, mut r) = fixture();
    r.transcript.length = usize::MAX;
    assert!(assess(&m, &r).is_err());
    r = request(&m, 19, 13);
    r.linearization.as_mut().unwrap().template_cut_offset = usize::MAX;
    assert!(assess(&m, &r).is_err());
    r = request(&m, 19, 13);
    r.transcript.start = usize::MAX;
    assert!(assess(&m, &r).is_err());
    r = request(&m, 19, 13);
    r.provenance.statement = "x".repeat(MAX_REQUEST_BYTES + 1);
    assert!(assess(&m, &r).is_err());
    let mut rna = m.clone();
    rna.sequence = Sequence::new(Alphabet::Rna, "AGGAAA").unwrap();
    r = request(&rna, 0, 6);
    assert!(assess(&rna, &r).is_err());
}
