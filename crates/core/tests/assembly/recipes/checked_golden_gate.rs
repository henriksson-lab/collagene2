//! Literal chemistry-aware cleavage and junction reassessment, independent of any target.
use collagene_core::{
    analysis,
    chemistry::{Duplex, Strand},
    golden_gate::{self, GoldenGateEnzyme, GoldenGateOptions},
    methylation::{MethylationState, Outcome},
    modified_duplex::{
        self as m, CallState, CheckedDigestOptions, Chemical, CleavageBudget,
        GenerationAssumptions, TrackedDuplex,
    },
};
use collagene_model::{Molecule, Topology};
fn assumptions(dcm: MethylationState) -> GenerationAssumptions {
    GenerationAssumptions {
        dam: MethylationState::Unmodified,
        dcm,
        cpg: MethylationState::Unmodified,
    }
}
fn donor(left: &str, payload: &str, right: &str) -> Molecule {
    Molecule::dna(
        "Literal donor",
        &format!("GGTCTCA{left}{payload}{right}AGAGACC"),
        Topology::Linear,
    )
    .unwrap()
}
fn options() -> GoldenGateOptions {
    GoldenGateOptions {
        max_products: 1000,
        max_fragments: 32,
        ..Default::default()
    }
}
fn generated(molecule: &Molecule, dcm: MethylationState) -> TrackedDuplex {
    TrackedDuplex::from_input(molecule, &assumptions(dcm)).unwrap()
}
#[test]
fn unknown_originals_withhold_the_entire_pool_but_actual_pcr_strands_support_literal_circle() {
    let sources = [donor("AATG", "CCCC", "GCTT"), donor("GCTT", "TTTT", "AATG")];
    let unknown = golden_gate::simulate(&sources, &options()).unwrap();
    assert!(!unknown.products_determined);
    assert!(!unknown.enumeration_complete);
    assert!(unknown.fragments.is_empty() && unknown.candidates.is_empty());
    assert_eq!(unknown.source_cleavage.len(), 2);
    assert!(
        unknown
            .source_cleavage
            .iter()
            .flat_map(|s| &s.sites)
            .any(|s| s.outcome == Outcome::Unknown)
    );
    let actual = sources
        .iter()
        .map(|s| TrackedDuplex::exact_pcr(&s.sequence).unwrap())
        .collect::<Vec<_>>();
    let mut budget = CleavageBudget::new(100_000_000).unwrap();
    let fresh =
        golden_gate::simulate_tracked(&sources, &actual, &options(), Some(0), &mut budget).unwrap();
    assert!(fresh.report.products_determined && fresh.report.enumeration_complete);
    assert_eq!(fresh.report.candidates.len(), 2);
    assert_eq!(fresh.report.comparisons, budget.used());
    assert!(
        fresh
            .report
            .candidates
            .iter()
            .any(|c| c.duplex.top.as_str() == "AATGCCCCGCTTTTTT")
    );
    assert_eq!(
        fresh.selected.unwrap().duplex,
        fresh.report.candidates[0].duplex
    );
    assert!(
        fresh
            .report
            .source_cleavage
            .iter()
            .flat_map(|s| &s.sites)
            .all(|s| s.outcome == Outcome::PermittedUnderProfile)
    );
}
#[test]
fn source_dcm_uncertainty_is_positional_not_a_generation_record_flag() {
    // Dcm CCAGG overlaps this BsaI GGTCTC at its first G (bottom methyl-C).
    // R0535's conflicting Dcm context is unsupported, never silently HFv2.
    let source = Molecule::dna("Dcm overlap", "AACCAGGTCTCAAAAAAAA", Topology::Circular).unwrap();
    let mut actual = generated(&source, MethylationState::FullyModified);
    assert!(m::summary(&actual).unwrap().known_modified > 0);
    actual.generation.clear(); // Native kernel unit test: chemical spans remain identical.
    let mut budget = CleavageBudget::new(100_000_000).unwrap();
    let result = golden_gate::simulate_tracked(
        std::slice::from_ref(&source),
        &[actual],
        &options(),
        None,
        &mut budget,
    )
    .unwrap()
    .report;
    assert!(
        !result.products_determined && result.fragments.is_empty() && result.candidates.is_empty()
    );
    assert_eq!(result.source_cleavage[0].sites[0].outcome, Outcome::Unknown);
    assert!(
        result.source_cleavage[0].sites[0]
            .evidence
            .iter()
            .any(|e| e.state == CallState::Modified)
    );
    assert_eq!(
        result.source_cleavage[0].options.selection.enzymes[0].id,
        "bsai"
    );
    assert_eq!(result.profile_revision, m::profile_revision());
}
#[test]
fn physically_inherited_methyl_c_at_a_new_junction_changes_recut_evidence_without_changing_sequence()
 {
    // Left retains top AATGCCAG and the first G of right's GTCT is supplied by
    // its cohesive boundary. The new BsaI site GGTCTC starts at output7.
    // Left's original CCAGG supplied a bottom methyl-C now inside that new site.
    let sources = [
        donor("AATG", "CCAG", "GTCT"),
        donor("GTCT", "CAAAA", "AATG"),
    ];
    let actual = sources
        .iter()
        .map(|s| generated(s, MethylationState::FullyModified))
        .collect::<Vec<_>>();
    let mut budget = CleavageBudget::new(100_000_000).unwrap();
    let report = golden_gate::simulate_tracked(&sources, &actual, &options(), None, &mut budget)
        .unwrap()
        .report;
    assert!(report.products_determined, "{:?}", report.source_cleavage);
    assert!(
        report
            .source_cleavage
            .iter()
            .flat_map(|s| &s.sites)
            .all(|s| s.outcome == Outcome::PermittedUnderProfile)
    );
    let index = report
        .candidates
        .iter()
        .position(|c| c.duplex.top.as_str() == "AATGCCAGGTCTCAAAA")
        .unwrap();
    let mut replay_budget = CleavageBudget::new(100_000_000).unwrap();
    let fresh = golden_gate::simulate_tracked(
        &sources,
        &actual,
        &options(),
        Some(index),
        &mut replay_budget,
    )
    .unwrap();
    let selected = fresh.selected.unwrap();
    let candidate = &fresh.report.candidates[index];
    assert_eq!(candidate.recut_sites.len(), 1);
    assert_eq!(candidate.recut_sites[0].recognition_start, 7);
    assert_eq!(candidate.recut_sites[0].outcome, Outcome::Unknown);
    assert!(!candidate.recut_sites[0].cleavable_under_model && !candidate.restriction_site_free);
    assert!(candidate.cleavage.complete && !candidate.cleavage.products_determined);
    assert!(
        candidate.cleavage.sites[0]
            .evidence
            .iter()
            .any(|e| e.chemical == Chemical::C5C && e.state == CallState::Modified)
    );
    assert_eq!(
        selected.state_at(Strand::Bottom, Chemical::C5C, 9).unwrap(),
        Some(CallState::Modified)
    );
    let independent_control = Molecule::dna(
        "Explicit unmodified circular control",
        "AATGCCAGGTCTCAAAA",
        Topology::Circular,
    )
    .unwrap();
    let copied = generated(&independent_control, MethylationState::Unmodified);
    assert_eq!(copied.duplex.top, selected.duplex.top);
    let mut copy_budget = CleavageBudget::new(100_000_000).unwrap();
    let (sites, evidence) =
        golden_gate::assess_product(&copied, GoldenGateEnzyme::BsaI, &mut copy_budget).unwrap();
    assert_eq!(sites[0].recognition_start, 7);
    assert_eq!(sites[0].outcome, Outcome::PermittedUnderProfile);
    assert!(evidence.products_determined && sites[0].cleavable_under_model);
    for reverse in [false, true] {
        let oriented = if reverse {
            selected.reverse().unwrap()
        } else {
            selected.clone()
        };
        let mut b = CleavageBudget::new(100_000_000).unwrap();
        let (sites, _) =
            golden_gate::assess_product(&oriented, GoldenGateEnzyme::BsaI, &mut b).unwrap();
        assert_eq!(sites.len(), 1);
        assert_eq!(sites[0].outcome, Outcome::Unknown);
    }
}
#[test]
fn blocked_intact_circle_survives_but_missing_exact_profile_never_borrows_by_geometry() {
    let stock =
        Molecule::dna("Actual full-CpG circle", "AAGCGGCCGCAA", Topology::Circular).unwrap();
    let blocked = TrackedDuplex::from_input(
        &stock,
        &GenerationAssumptions {
            dam: MethylationState::Unmodified,
            dcm: MethylationState::Unmodified,
            cpg: MethylationState::FullyModified,
        },
    )
    .unwrap();
    let o = CheckedDigestOptions {
        selection: m::checked_selection("NotI").unwrap(),
        ..Default::default()
    };
    let r = m::checked_digest(&blocked, &o).unwrap();
    assert!(r.complete && r.products_determined);
    assert_eq!(r.sites[0].outcome, Outcome::Blocked);
    assert_eq!(r.fragments.len(), 1);
    assert_eq!(r.fragments[0].duplex.duplex, blocked.duplex);
    // AatII isn't a starter positional profile, but actual unmodified strands
    // use its pinned geometry rather than an unrelated sensitivity assumption.
    assert!(
        !analysis::starter_enzymes()
            .iter()
            .any(|e| e.name == "AatII")
    );
    let stock = Molecule::dna("AatII literal10", "AAGACGTCAA", Topology::Linear).unwrap();
    let clean = TrackedDuplex::exact_pcr(&stock.sequence).unwrap();
    let mut o = CheckedDigestOptions {
        selection: m::checked_selection("AatII").unwrap(),
        ..Default::default()
    };
    let r = m::checked_digest(&clean, &o).unwrap();
    assert!(r.products_determined);
    assert_eq!(r.fragments[0].duplex.duplex.top.as_str(), "AAGACGT");
    assert_eq!(r.fragments[0].duplex.duplex.bottom.as_str(), "CTT");
    let unknown =
        TrackedDuplex::from_unknown_duplex(&Duplex::from_molecule(&stock).unwrap()).unwrap();
    let r = m::checked_digest(&unknown, &o).unwrap();
    assert!(!r.products_determined && r.fragments.is_empty());
    o.selection.enzymes[0].enzyme.top_cut += 1;
    let mut budget = CleavageBudget::new(100_000_000).unwrap();
    assert!(m::checked_digest_with_budget(&clean, &o, &mut budget).is_err());
    assert_eq!(budget.used(), 0);
}
#[test]
fn preparation_source_and_recut_scans_share_one_ledger_and_capped_recut_retains_actual_candidate() {
    let sources = [donor("AATG", "CCCC", "GCTT"), donor("GCTT", "TTTT", "AATG")];
    let actual = sources
        .iter()
        .map(|s| TrackedDuplex::exact_pcr(&s.sequence).unwrap())
        .collect::<Vec<_>>();
    let mut unlimited = CleavageBudget::new(100_000_000).unwrap();
    let full = golden_gate::simulate_tracked(&sources, &actual, &options(), None, &mut unlimited)
        .unwrap()
        .report;
    let source_work = full
        .source_cleavage
        .iter()
        .map(|s| s.comparisons)
        .sum::<usize>();
    assert!(full.comparisons > source_work);
    let mut limited = CleavageBudget::new(source_work).unwrap();
    let capped =
        golden_gate::simulate_tracked(&sources, &actual, &options(), Some(0), &mut limited)
            .unwrap();
    assert!(capped.report.products_determined && !capped.report.enumeration_complete);
    assert!(capped.selected.is_some());
    assert!(
        capped
            .report
            .candidates
            .iter()
            .all(|c| !c.restriction_site_free && !c.cleavage.complete)
    );
    assert_eq!(capped.report.comparisons, source_work);
    let o = CheckedDigestOptions {
        selection: m::checked_selection("BsaI").unwrap(),
        ..Default::default()
    };
    let prep = m::checked_digest(&actual[0], &o).unwrap();
    let mut limited = CleavageBudget::new(prep.comparisons + 1).unwrap();
    let first = m::checked_digest_with_budget(&actual[0], &o, &mut limited).unwrap();
    assert!(first.products_determined);
    let withheld = golden_gate::simulate_tracked(&sources, &actual, &options(), None, &mut limited)
        .unwrap()
        .report;
    assert!(
        !withheld.products_determined
            && withheld.candidates.is_empty()
            && withheld.fragments.is_empty()
    );
    assert_eq!(withheld.comparisons, limited.used());
    assert_eq!(limited.remaining(), 0);
}

#[test]
fn exact_recognition_budget_proves_absence_without_irrelevant_context_work() {
    use collagene_core::restriction_inspection::{self, InspectionOptions};
    let source = Molecule::dna(
        "No selected motif, unknown chemical state",
        "AACCGGTT",
        Topology::Circular,
    )
    .unwrap();
    let actual =
        TrackedDuplex::from_unknown_duplex(&Duplex::from_molecule(&source).unwrap()).unwrap();
    let scan = restriction_inspection::scan(
        "control",
        &source,
        &[GoldenGateEnzyme::BsaI.profile()],
        &InspectionOptions {
            max_comparisons: 100_000_000,
            max_sites: 4096,
            max_fragments: 4097,
            max_report_bytes: 8 * 1024 * 1024,
        },
    )
    .unwrap();
    assert!(scan.scan_complete && scan.sites.is_empty());
    let exact = scan.comparisons as usize;
    let mut budget = CleavageBudget::new(exact).unwrap();
    let fresh = golden_gate::simulate_tracked(
        std::slice::from_ref(&source),
        std::slice::from_ref(&actual),
        &options(),
        None,
        &mut budget,
    )
    .unwrap();
    assert_eq!(budget.remaining(), 0);
    assert_eq!(fresh.report.comparisons, exact);
    assert!(fresh.report.products_determined && fresh.report.enumeration_complete);
    assert_eq!(fresh.report.carryover_sources, [0]);
    assert!(fresh.report.candidates.is_empty());
    assert!(fresh.report.source_cleavage[0].complete);
    assert!(fresh.report.source_cleavage[0].limits.is_empty());
    let mut budget = CleavageBudget::new(exact).unwrap();
    let (sites, recut) =
        golden_gate::assess_product(&actual, GoldenGateEnzyme::BsaI, &mut budget).unwrap();
    assert!(sites.is_empty() && recut.complete && recut.products_determined);
    assert_eq!(recut.comparisons, exact);
    assert_eq!(budget.remaining(), 0);
}

#[test]
fn actual_unmodified_recognition_residues_need_no_context_budget() {
    use collagene_core::restriction_inspection::{self, InspectionOptions};
    let source = Molecule::dna(
        "Known new PCR with BsaI site",
        "AAGGTCTCAAAAACCCCTTTT",
        Topology::Linear,
    )
    .unwrap();
    let actual = TrackedDuplex::exact_pcr(&source.sequence).unwrap();
    let scan = restriction_inspection::scan(
        "control",
        &source,
        &[GoldenGateEnzyme::BsaI.profile()],
        &InspectionOptions {
            max_comparisons: 100_000_000,
            max_sites: 4096,
            max_fragments: 4097,
            max_report_bytes: 8 * 1024 * 1024,
        },
    )
    .unwrap();
    assert!(scan.scan_complete);
    assert_eq!(scan.sites.len(), 1);
    let mut budget = CleavageBudget::new(scan.comparisons as usize).unwrap();
    let report = m::checked_digest_with_budget(
        &actual,
        &CheckedDigestOptions {
            selection: m::checked_selection("BsaI").unwrap(),
            max_comparisons: scan.comparisons as usize,
            max_sites: 4096,
        },
        &mut budget,
    )
    .unwrap();
    assert!(report.complete && report.products_determined);
    assert_eq!(report.fragments.len(), 2);
    assert_eq!(report.sites[0].outcome, Outcome::PermittedUnderProfile);
    assert_eq!(budget.remaining(), 0);
    assert_eq!(report.comparisons, scan.comparisons as usize);
}
