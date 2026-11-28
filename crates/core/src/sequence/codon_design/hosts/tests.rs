use super::*;
use std::io::Write;

fn metadata() -> Dataset {
    let mut d = Dataset::mg1655();
    d.id = "synthetic-derivation-fixture".into();
    d
}
fn report() -> String {
    format!(
        "# Assembly name: ASM584v2\n# Organism name: Escherichia coli str. K-12 substr. MG1655\n# Taxid: 511145\n# RefSeq assembly accession: GCF_000005845.2\nU00096.3\tassembled-molecule\tna\tChromosome\tU00096.3\t=\tNC_000913.3\tPrimary Assembly\t4641652\tna\n"
    )
}
fn annotation() -> &'static [u8] {
    b"# Assembly accession\tFeatures changed\nGCF_000005845.2\tfixture\n"
}
fn gff_prefix() -> &'static str {
    "##gff-version 3\n#!genome-build-accession NCBI_Assembly:GCF_000005845.2\n##sequence-region NC_000913.3 1 4641652\n"
}
fn row(locus: &str, protein: &str, start: usize, bases: &str, flags: &str) -> (String, String) {
    let end = start + bases.len() - 1;
    (
        format!(
            ">lcl|NC_000913.3_cds_{protein}_1 [locus_tag={locus}] [protein=fixture [4Fe-4S] protein] [protein_id={protein}] [location={start}..{end}] [gbkey=CDS]\n{bases}\n"
        ),
        format!(
            "NC_000913.3\tRefSeq\tgene\t{start}\t{end}\t.\t+\t.\tID=gene-{locus};locus_tag={locus}\nNC_000913.3\tRefSeq\tCDS\t{start}\t{end}\t.\t+\t0\tID=cds-{protein}-{locus};Parent=gene-{locus};locus_tag={locus};protein_id={protein};transl_table=11{flags}\n"
        ),
    )
}
fn derive(fasta: &str, gff: &str) -> Result<HostProfile, Error> {
    HostProfile::derive_expanded(
        metadata(),
        Expanded {
            cds: fasta.as_bytes(),
            gff: gff.as_bytes(),
            assembly_report: report().as_bytes(),
            annotation_hashes: annotation(),
        },
    )
}
#[test]
fn independently_specified_internal_counts_start_stop_and_nonzero_weights() {
    let (f1, g1) = row("locus-a", "protein-a", 1, "ATGGCTGCTTAA", "");
    let (f2, g2) = row("locus-b", "protein-b", 101, "GTGGCCGCTTAG", "");
    let profile = derive(&(f1 + &f2), &(gff_prefix().to_owned() + &g1 + &g2)).unwrap();
    assert_eq!(profile.observed_loci, 2);
    assert_eq!(profile.admitted_loci, 2);
    assert_eq!(profile.internal_codons, 4);
    assert_eq!(profile.counts[index(b"GCT").unwrap()], 3);
    assert_eq!(profile.counts[index(b"GCC").unwrap()], 1);
    assert_eq!(profile.initiation_counts[index(b"ATG").unwrap()], 1);
    assert_eq!(profile.initiation_counts[index(b"GTG").unwrap()], 1);
    assert_eq!(profile.terminal_stop_counts[index(b"TAA").unwrap()], 1);
    assert_eq!(profile.terminal_stop_counts[index(b"TAG").unwrap()], 1);
    let usage = profile.usage().unwrap();
    assert_eq!(
        usage[index(b"GCT").unwrap()].synonymous_millionths,
        Some(750_000)
    );
    assert_eq!(usage[index(b"GCC").unwrap()].preference, 333_334); // 1 + floor(999999 / 3)
    assert_eq!(usage[index(b"GCT").unwrap()].preference, 1_000_000);
    assert_eq!(usage[index(b"GCA").unwrap()].preference, 1);
    assert!(profile.preferences().unwrap().iter().all(|w| w.weight > 0));
    for stop in [b"TAA", b"TAG", b"TGA"] {
        assert_eq!(usage[index(stop).unwrap()].preference, 1);
    }
}
#[test]
fn true_paralog_loci_are_counted_without_protein_or_sequence_deduplication() {
    let (f1, g1) = row("copy-a", "same-protein", 1, "ATGAAATAA", "");
    let (f2, g2) = row("copy-b", "same-protein", 101, "ATGAAATAA", "");
    let a = derive(&(f1.clone() + &f2), &(gff_prefix().to_owned() + &g1 + &g2)).unwrap();
    let b = derive(&(f2 + &f1), &(gff_prefix().to_owned() + &g2 + &g1)).unwrap();
    assert_eq!(a.counts[index(b"AAA").unwrap()], 2);
    assert_eq!(a.admitted_loci, 2);
    assert_eq!(a.counts, b.counts);
    assert_eq!(a.admitted_loci_blake3, b.admitted_loci_blake3);
    assert_ne!(
        a.identity, b.identity,
        "raw byte order belongs to source provenance"
    );
}
#[test]
fn scientific_exclusions_are_mutually_exclusive_and_accounted() {
    let cases = [
        ("valid", "ATGAAATAA", ""),
        ("pseudo", "ATGAAATAA", ";pseudo=true"),
        ("exception", "ATGAAATAA", ";transl_except=fixture"),
        ("partial", "ATGAAATAA", ";partial=true"),
        ("ambiguous", "ATGNNNTAA", ""),
        ("start", "CCCAAATAA", ""),
        ("terminal", "ATGAAAAAA", ""),
        ("internal", "ATGTAATAA", ""),
    ];
    let mut fasta = String::new();
    let mut gff = gff_prefix().to_owned();
    for (i, (locus, bases, flags)) in cases.iter().enumerate() {
        let (f, g) = row(locus, &format!("protein-{i}"), 1 + i * 100, bases, flags);
        fasta += &f;
        gff += &g;
    }
    let profile = derive(&fasta, &gff).unwrap();
    assert_eq!(profile.admitted_loci, 1);
    assert_eq!(profile.exclusions.iter().map(|r| r.loci).sum::<u64>(), 7);
    assert_eq!(profile.exclusions.len(), 7);
}
#[test]
fn mismatched_duplicate_or_missing_locus_coordinates_and_codes_fail() {
    let (f, g) = row("locus", "protein", 1, "ATGAAATAA", "");
    let g = gff_prefix().to_owned() + &g;
    for invalid in [
        f.replace("locus_tag=locus", "locus_tag=other"),
        f.replace("location=1..9", "location=2..10"),
        f.clone() + &f,
        f.replace("[gbkey=CDS]", "[transl_table=1] [gbkey=CDS]"),
        f.replace("ATGAAATAA", "ATGAAATAZ"),
    ] {
        assert!(derive(&invalid, &g).is_err());
    }
    assert!(derive(&f, &g.replace("GCF_000005845.2", "GCF_000005845.3")).is_err());
    assert!(derive(&f, &(g.clone() + g.lines().last().unwrap() + "\n")).is_err());
    assert!(derive("", &g).is_err());
}
#[test]
fn malformed_and_resource_limits_never_produce_partial_profile() {
    let (f, g) = row("locus", "protein", 1, "ATGAAATAA", "");
    let g = gff_prefix().to_owned() + &g;
    let huge = f.replace("ATGAAATAA", &"A".repeat(MAX_CDS_BASES + 1));
    assert!(matches!(derive(&huge, &g), Err(Error::TooLarge)));
    assert!(derive(&f, &format!("{g}{}\n", "#".repeat(MAX_LINE_BYTES + 1))).is_err());
    let gzip = |bytes: &[u8]| {
        let mut e = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
        e.write_all(bytes).unwrap();
        e.finish().unwrap()
    };
    let fz = gzip(f.as_bytes());
    let gz = gzip(g.as_bytes());
    let report = report();
    let profile = HostProfile::derive_gzip(
        metadata(),
        Compressed {
            cds_gzip: &fz,
            gff_gzip: &gz,
            assembly_report: report.as_bytes(),
            annotation_hashes: annotation(),
        },
    )
    .unwrap();
    assert_eq!(profile.internal_codons, 1);
    let mut truncated = fz.clone();
    truncated.pop();
    assert!(
        HostProfile::derive_gzip(
            metadata(),
            Compressed {
                cds_gzip: &truncated,
                gff_gzip: &gz,
                assembly_report: report.as_bytes(),
                annotation_hashes: annotation()
            }
        )
        .is_err()
    );
    let mut junk = fz;
    junk.push(42);
    assert!(
        HostProfile::derive_gzip(
            metadata(),
            Compressed {
                cds_gzip: &junk,
                gff_gzip: &gz,
                assembly_report: report.as_bytes(),
                annotation_hashes: annotation()
            }
        )
        .is_err()
    );
}
#[test]
fn profile_identity_and_conservation_checks_reject_forgery() {
    let (f, g) = row("locus", "protein", 1, "ATGAAATAA", "");
    let profile = derive(&f, &(gff_prefix().to_owned() + &g)).unwrap();
    let mut forged = profile.clone();
    forged.counts[index(b"AAA").unwrap()] += 1;
    assert!(forged.validate().is_err());
    let mut forged = profile.clone();
    forged.dataset.assembly_name = "other".into();
    assert!(forged.validate().is_err());
    let mut forged = profile.clone();
    forged.version += 1;
    assert!(matches!(forged.validate(), Err(Error::Unsupported(_))));
    let mut forged = profile;
    forged.counts[index(b"AAA").unwrap()] = u64::MAX;
    assert!(forged.validate().is_err());
}
#[test]
fn captured_mg1655_cds_gff_release_derives_freshly_with_complete_accounting() {
    let profile = HostProfile::derive_mg1655().unwrap();
    profile.validate().unwrap();
    assert_eq!(profile.observed_loci, 4318);
    assert_eq!(
        profile
            .exclusions
            .iter()
            .find(|r| r.reason == Exclusion::Pseudogene)
            .unwrap()
            .loci,
        18
    );
    assert_eq!(
        profile
            .exclusions
            .iter()
            .find(|r| r.reason == Exclusion::TranslationException)
            .unwrap()
            .loci,
        6
    );
    assert_eq!(
        profile.admitted_loci + profile.exclusions.iter().map(|r| r.loci).sum::<u64>(),
        4318
    );
    assert!(profile.annotation_hashes.contains("2026-09-04 17:13:00"));
    assert_eq!(profile.preferences().unwrap().len(), 64);
    assert!(profile.sources[0].compressed);
    assert_eq!(profile.sources[0].expanded_bytes, 5_043_083);
    // Capture with --nocapture for the parent integration's one-time catalogue
    // emission; this is real freshly derived output, never placeholder counts.
    println!(
        "MG1655_PROFILE_JSON={}",
        serde_json::to_string(&profile).unwrap()
    );
}

#[test]
fn compact_runtime_catalogue_matches_all_independent_counts_and_fresh_derivation() {
    let expected: serde_json::Value =
        serde_json::from_str(include_str!("data/mg1655/expected-counts.json")).unwrap();
    let summary = catalogue().unwrap().remove(0);
    let profile = resolve(&summary.reference).unwrap();
    assert_eq!(profile, HostProfile::derive_mg1655().unwrap());
    assert_eq!(
        serde_json::to_value(&profile.counts).unwrap(),
        expected["counts"]
    );
    assert_eq!(
        serde_json::to_value(&profile.initiation_counts).unwrap(),
        expected["initiation_counts"]
    );
    assert_eq!(
        serde_json::to_value(&profile.terminal_stop_counts).unwrap(),
        expected["terminal_stop_counts"]
    );
    assert_eq!(summary.admitted_loci, 4294);
    assert_eq!(summary.excluded_loci, 24);
    assert_eq!(summary.internal_codons, 1322130);
    let mut wrong = summary.reference;
    wrong.id.push_str("-unavailable");
    assert!(resolve(&wrong).is_err());
}
#[test]
fn pinned_preferences_are_exclusive_nonzero_bound_and_retained_as_documentary_evidence() {
    use crate::sequence::codon_design::{self, CodonPreferences, DesignOptions};
    use collagene_model::{Alphabet, Molecule, Sequence, Topology};
    let summary = catalogue().unwrap().remove(0);
    let source = Molecule {
        name: "Required peptide".into(),
        sequence: Sequence::new(Alphabet::Protein, "MK*").unwrap(),
        topology: Topology::Linear,
        features: vec![],
        metadata: Default::default(),
        original_genbank: None,
    };
    let options = DesignOptions {
        genetic_code: 11,
        preferences: CodonPreferences::ReferenceHost {
            reference: summary.reference.clone(),
        },
        ..Default::default()
    };
    let profile = resolve(&summary.reference).unwrap();
    assert!(profile.preferences().unwrap().iter().all(|w| w.weight > 0));
    let report = codon_design::design(&source, &options).unwrap();
    report.validate().unwrap();
    let candidate = report.candidate.as_ref().unwrap();
    assert_eq!(candidate.sequence.as_str(), "ATGAAATAA");
    assert_eq!(report.score, Some(2_000_001));
    let checked = codon_design::check(&source, candidate, &options).unwrap();
    checked.validate_for(&source, candidate, &options).unwrap();
    assert_eq!(
        checked.preference_evidence.reference_host.as_ref(),
        Some(&profile)
    );
    assert!(candidate.metadata["codon_design:preference_evidence"].contains(&profile.identity));
    let explicit = DesignOptions {
        preferences: CodonPreferences::Explicit {
            weights: profile.preferences().unwrap(),
        },
        ..options.clone()
    };
    assert_ne!(
        codon_design::options_revision(&options).unwrap(),
        codon_design::options_revision(&explicit).unwrap()
    );
    assert_eq!(
        codon_design::design(&source, &explicit).unwrap().score,
        report.score
    );
    let mut forged = report.clone();
    forged.preference_evidence.expanded_weights[0].weight = 0;
    assert!(forged.validate().is_err());
    let mut incompatible = options.clone();
    incompatible.genetic_code = 1;
    assert!(incompatible.validate().is_err());
    let mut missing = options;
    missing.preferences = CodonPreferences::Explicit { weights: vec![] };
    assert!(missing.validate().is_err());
}
