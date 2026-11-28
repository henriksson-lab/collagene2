use collagene_core::{
    formats::sequence::*,
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    protocol::{
        execution::process_request,
        load_document,
        request::Request,
        save_document,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::analysis::*,
    verification::{engine::verify, evidence::Status, output::export_protocol_html},
};
use collagene_model::{Alphabet, Feature, Location, Molecule, Sequence, Topology};

fn molecule(bases: &str, topology: Topology) -> Molecule {
    Molecule::dna("test", bases, topology).unwrap()
}
fn feature(location: Location) -> Feature {
    Feature {
        kind: "CDS".into(),
        location,
        qualifiers: vec![
            ("label".into(), Some("test gene".into())),
            ("translation".into(), Some("stale".into())),
        ],
        needs_review: false,
    }
}

#[test]
fn alphabets_normalize_without_inventing_bases() {
    assert_eq!(Sequence::dna(" aCgT\nN ").unwrap().as_str(), "ACGTN");
    for bad in ["AC?T", "ACUT", "12AC", "AβC", "---"] {
        assert!(Sequence::dna(bad).is_err(), "{bad}");
    }
    assert!(Sequence::dna("  ").is_err());
    assert_eq!(
        Sequence::new(Alphabet::Rna, "augc")
            .unwrap()
            .reverse_complement()
            .unwrap()
            .as_str(),
        "GCAU"
    );
    assert!(
        Sequence::new(Alphabet::Protein, "MEK*")
            .unwrap()
            .reverse_complement()
            .is_err()
    );
}

#[test]
fn iupac_complement_is_an_involution() {
    let dna = Sequence::dna("ACGTRYSWKMBDHVN").unwrap();
    assert_eq!(
        dna.reverse_complement().unwrap().as_str(),
        "NBDHVKMWSRYACGT"
    );
    assert_eq!(
        dna.reverse_complement()
            .unwrap()
            .reverse_complement()
            .unwrap(),
        dna
    );
    assert_eq!(Sequence::dna("GCNNAT").unwrap().gc_fraction(), Some(0.5));
    assert_eq!(Sequence::dna("NNN").unwrap().gc_fraction(), None);
}

#[test]
fn exact_compound_and_reverse_location_extraction() {
    let sequence = Sequence::dna("AAAACCCCGGGGTTTT").unwrap();
    let origin = Location::Join(vec![Location::range(12, 16), Location::range(0, 4)]);
    assert_eq!(origin.extract(&sequence).unwrap(), "TTTTAAAA");
    assert_eq!(
        Location::Complement(Box::new(Location::Join(vec![
            Location::range(0, 4),
            Location::range(4, 8)
        ])))
        .extract(&sequence)
        .unwrap(),
        "GGGGTTTT"
    );
    assert!(
        Location::Order(vec![Location::range(0, 4)])
            .extract(&sequence)
            .is_err()
    );
    assert!(
        Location::Range {
            start: 0,
            end: 4,
            before: true,
            after: false
        }
        .extract(&sequence)
        .is_err()
    );
}

#[test]
fn edit_remaps_coordinates_and_invalidates_translation() {
    let mut mol = molecule("AAAACCCCGGGG", Topology::Circular);
    mol.features = vec![
        feature(Location::range(4, 8)),
        feature(Location::Join(vec![
            Location::range(8, 12),
            Location::range(0, 2),
        ])),
    ];
    let edited = mol.replace(5, 7, "TTTT").unwrap();
    assert_eq!(edited.sequence.as_str(), "AAAACTTTTCGGGG");
    assert_eq!(edited.features[0].location, Location::range(4, 10));
    assert!(edited.features[0].needs_review);
    assert!(
        !edited.features[0]
            .qualifiers
            .iter()
            .any(|(k, _)| k == "translation")
    );
    assert_eq!(
        edited.features[1].location,
        Location::Join(vec![Location::range(10, 14), Location::range(0, 2)])
    );
    assert!(!edited.features[1].needs_review);
    assert_eq!(mol.sequence.as_str(), "AAAACCCCGGGG");
    assert_ne!(mol.revision(), edited.revision());
}

#[test]
fn edit_boundary_policy_does_not_invent_feature_extensions() {
    let mut mol = molecule("AAAACCCC", Topology::Linear);
    mol.features.push(feature(Location::range(4, 8)));
    assert_eq!(
        mol.replace(4, 4, "GG").unwrap().features[0].location,
        Location::range(6, 10)
    );
    assert_eq!(
        mol.replace(8, 8, "GG").unwrap().features[0].location,
        Location::range(4, 8)
    );
    assert!(mol.replace(0, 8, "").is_err());
    assert!(mol.replace(4, 8, "").unwrap().features.is_empty());
}

#[test]
fn unsupported_coordinate_metadata_is_not_exported_as_if_remapped() {
    let mut mol = molecule("AAAACCCC", Topology::Linear);
    let mut f = feature(Location::range(0, 8));
    f.qualifiers
        .push(("transl_except".into(), Some("(pos:4..6,aa:Sec)".into())));
    mol.features.push(f);
    let original = mol.clone();
    let edited = mol.replace(0, 0, "AA").unwrap();
    assert_eq!(edited.sequence.as_str(), "AAAAAACCCC");
    assert_eq!(edited.features[0].location, Location::range(2, 10));
    assert!(!edited.features[0].needs_review);
    assert!(
        edited.features[0]
            .qualifiers
            .contains(&("transl_except".into(), Some("(pos:6..8,aa:Sec)".into())))
    );
    assert_eq!(mol, original);
    mol.features[0]
        .qualifiers
        .push(("rpt_unit_range".into(), Some("4..6".into())));
    assert!(mol.replace(0, 0, "AA").is_err());
    let mut mol = molecule("AAAACCCC", Topology::Linear);
    mol.original_genbank = Some("REFERENCE   1 (bases 1 to 8)\n".into());
    assert!(mol.replace(0, 0, "AA").is_err());
}

#[test]
fn circular_search_sees_origin_and_flags_ambiguity() {
    let mol = molecule("ATTCCCGA", Topology::Circular);
    let hits = find_matches(&mol, "GAATTC", false).unwrap();
    assert_eq!(
        hits,
        vec![Match {
            start: 6,
            length: 6,
            reverse: false,
            certain: true
        }]
    );
    let mol = molecule("GNATTC", Topology::Linear);
    assert!(!find_matches(&mol, "GAATTC", true).unwrap()[0].certain);
    assert!(
        find_matches(&molecule("AGCT", Topology::Linear), "AGCT", true)
            .unwrap()
            .len()
            == 1
    );
}

#[test]
fn restriction_cut_offsets_and_fragment_accounting() {
    let eco = starter_enzymes().remove(0);
    let mol = molecule("CCGAATTCAA", Topology::Linear);
    let sites = restriction_sites(&mol, &eco).unwrap();
    assert_eq!((sites[0].top_cut, sites[0].bottom_cut), (3, 7));
    assert_eq!(digest_lengths(&mol, &sites).unwrap(), vec![7, 3]);
    let circle = molecule("ATTCCCGA", Topology::Circular);
    let sites = restriction_sites(&circle, &eco).unwrap();
    assert_eq!((sites[0].top_cut, sites[0].bottom_cut), (7, 3));
    assert_eq!(digest_lengths(&circle, &sites).unwrap(), vec![8]);
    let mut two = sites.clone();
    let mut second = sites[0].clone();
    second.top_cut = 2;
    two.push(second);
    assert_eq!(digest_lengths(&circle, &two).unwrap(), vec![5, 3]);
}

#[test]
fn type_iis_cuts_are_mirrored_on_reverse_strand() {
    let bsa = starter_enzymes()
        .into_iter()
        .find(|e| e.name == "BsaI")
        .unwrap();
    let forward = molecule("AAGGTCTCAAAACCCCGGGG", Topology::Linear);
    let sites = restriction_sites(&forward, &bsa).unwrap();
    assert_eq!((sites[0].top_cut, sites[0].bottom_cut), (9, 13));
    let reverse = molecule(
        forward.sequence.reverse_complement().unwrap().as_str(),
        Topology::Linear,
    );
    let rev = restriction_sites(&reverse, &bsa).unwrap();
    assert!(rev[0].reverse);
    assert_eq!(
        (rev[0].top_cut, rev[0].bottom_cut),
        (forward.sequence.len() - 13, forward.sequence.len() - 9)
    );
}

#[test]
fn type_iis_sites_without_complete_cut_context_are_not_silently_omitted() {
    let bsa = starter_enzymes()
        .into_iter()
        .find(|e| e.name == "BsaI")
        .unwrap();
    assert!(restriction_sites(&molecule("GGTCTCAA", Topology::Linear), &bsa).is_err());
}

#[test]
fn exact_pcr_includes_both_tails_in_correct_orientation() {
    let mol = molecule("ACGACCCCGGGGTAAC", Topology::Linear);
    let f = Primer::new("F", "ACGA", "CGT").unwrap();
    let r = Primer::new("R", "GTTA", "GAC").unwrap();
    let products = simulate_pcr(&mol, &f, &r).unwrap();
    assert_eq!(products.len(), 1);
    assert_eq!(products[0].sequence.as_str(), "CGTACGACCCCGGGGTAACGTC");
    assert!(!products[0].crosses_origin);
}

#[test]
fn pcr_crosses_circular_origin_and_rejects_unknown_binding() {
    let circle = molecule("CCCCTAACGGGGACGA", Topology::Circular);
    let f = Primer::new("F", "ACGA", "").unwrap();
    let r = Primer::new("R", "GTTA", "").unwrap();
    let products = simulate_pcr(&circle, &f, &r).unwrap();
    assert_eq!(products.len(), 1);
    assert!(products[0].crosses_origin);
    assert_eq!(products[0].sequence.as_str(), "ACGACCCCTAAC");
    let unknown = molecule("ACGACCCCGGGGTAAN", Topology::Linear);
    assert!(simulate_pcr(&unknown, &f, &r).is_err());
}

#[test]
fn pcr_enumerates_competing_and_single_oligo_products() {
    let mol = molecule("AAAACCCCGGGGTTTT", Topology::Linear);
    let f = Primer::new("F", "AAAA", "CGT").unwrap();
    let r = Primer::new("R", "AAAA", "GAC").unwrap();
    let products = simulate_pcr(&mol, &f, &r).unwrap();
    assert_eq!(products.len(), 4);
    assert!(
        products
            .iter()
            .any(|p| p.forward_primer == "F" && p.reverse_primer == "F")
    );
    assert!(
        products
            .iter()
            .any(|p| p.forward_primer == "R" && p.reverse_primer == "F")
    );
}

#[test]
fn rotation_comparison_matches_naive_oracle() {
    for bits in 0..256u32 {
        let s: String = (0..8)
            .map(|i| if bits & (1 << i) == 0 { 'A' } else { 'C' })
            .collect();
        for shift in 0..8 {
            let b = format!("{}{}", &s[shift..], &s[..shift]);
            let expected = (0..8).find(|&i| format!("{}{}", &s[i..], &s[..i]) == b);
            assert_eq!(rotation_offset(&s, &b), expected);
        }
    }
    assert_eq!(rotation_offset("ACGT", "ACGA"), None);
    assert_eq!(rotation_offset("ACGT", "ACG"), None);
}

const GB: &str = "LOCUS       fixture                   16 bp    DNA     circular SYN 01-JAN-2026\nDEFINITION  Independent synthetic coordinate fixture.\nFEATURES             Location/Qualifiers\n     CDS             complement(join(13..16,1..4))\n                     /label=\"origin feature\"\n                     /note=\"first\"\n                     /note=\"second\"\n                     /ApEinfo_fwdcolor=\"#00aa77\"\nORIGIN\n        1 aaaaccccgg ggtttt\n//\n";

#[test]
fn genbank_round_trip_preserves_compound_locations_and_qualifiers() {
    let mol = import_genbank(GB).unwrap().molecules.remove(0);
    assert_eq!(mol.topology, Topology::Circular);
    assert_eq!(
        mol.features[0].location.extract(&mol.sequence).unwrap(),
        "TTTTAAAA"
    );
    let exported = export_genbank(&mol).unwrap();
    let round = import_genbank(&exported).unwrap().molecules.remove(0);
    assert_eq!(mol.features, round.features);
    assert_eq!(mol.sequence, round.sequence);
    assert_eq!(
        round.features[0]
            .qualifiers
            .iter()
            .filter(|(k, _)| k == "note")
            .count(),
        2
    );
}

#[test]
fn genbank_location_recursion_is_bounded_before_parsing() {
    let deep = format!("{}1..16{}", "complement(".repeat(40), ")".repeat(40));
    let text = GB.replace("complement(join(13..16,1..4))", &deep);
    assert!(
        import_genbank(&text)
            .unwrap_err()
            .to_string()
            .contains("nesting")
    );
}

#[test]
fn multi_record_genbank_retains_each_source_once() {
    let text = format!("{GB}{}", GB.replace("fixture", "second"));
    let report = import_genbank(&text).unwrap();
    assert_eq!(report.molecules.len(), 2);
    let retained: usize = report
        .molecules
        .iter()
        .map(|m| m.original_genbank.as_ref().unwrap().len())
        .sum();
    assert_eq!(retained, text.len());
    for molecule in &report.molecules {
        let output = export_genbank(molecule).unwrap();
        let round = import_genbank(&output).unwrap();
        assert_eq!(round.molecules[0].name, molecule.name);
    }
}

#[test]
fn fasta_rejects_empty_and_malformed_records_and_exports_losslessly() {
    assert!(import_fasta(">empty\n>next\nACGT").is_err());
    assert!(import_fasta("ACGT\n>x\nACGT").is_err());
    let docs = import_fasta(">first\nacgt\n>second\nNNNN").unwrap();
    assert_eq!(docs.molecules.len(), 2);
    let mol = &docs.molecules[0];
    let round = import_fasta(&export_fasta(mol)).unwrap();
    assert_eq!(round.molecules[0].sequence, mol.sequence);
}

#[test]
fn protein_export_uses_explicit_fasta_alphabet_instead_of_incorrect_genbank_units() {
    let mut mol = molecule("ACGT", Topology::Linear);
    mol.sequence = Sequence::new(Alphabet::Protein, "MEK*").unwrap();
    assert!(export_genbank(&mol).is_err());
    let fasta = export_fasta(&mol);
    let round = import_fasta(&fasta).unwrap().molecules.remove(0);
    assert_eq!(round.sequence, mol.sequence);
    assert_eq!(
        export_fasta(&round).matches("[molecule=protein]").count(),
        1
    );
}

#[test]
fn native_file_checksums_versions_and_semantics() {
    let mol = import_genbank(GB).unwrap().molecules.remove(0);
    let bytes = save_document(&mol).unwrap();
    assert_eq!(load_document(&bytes).unwrap(), mol);
    for i in [0, 4, 6, 10, bytes.len() - 1] {
        let mut corrupt = bytes.clone();
        corrupt[i] ^= 1;
        assert!(load_document(&corrupt).is_err());
    }
    assert!(load_document(&bytes[..bytes.len() - 1]).is_err());
    let mut invalid = mol;
    invalid.features[0].location = Location::range(0, 999);
    assert!(save_document(&invalid).is_err());
}

#[test]
fn location_decoding_limits_nesting_and_recovers_after_an_error() {
    let mut location = Location::range(0, 1);
    for _ in 0..40 {
        location = Location::Complement(Box::new(location));
    }
    let bytes = postcard::to_allocvec(&location).unwrap();
    assert!(postcard::from_bytes::<Location>(&bytes).is_err());
    let simple = Location::range(0, 1);
    let bytes = postcard::to_allocvec(&simple).unwrap();
    assert_eq!(postcard::from_bytes::<Location>(&bytes).unwrap(), simple);
}

#[test]
fn binary_worker_messages_bind_results_to_inputs() {
    let mol = molecule("AAAAGAATTC", Topology::Linear);
    let request = Request {
        id: 7,
        revision: mol.revision(),
        job: AnalysisJob::Search {
            molecule: mol,
            query: "GAATTC".into(),
        },
    };
    let response = collagene_core::protocol::response::Response::decode(
        &process_request(&collagene_core::protocol::request::Request::encode(&request).unwrap())
            .unwrap(),
    )
    .unwrap();
    assert_eq!(response.id, 7);
    assert_eq!(response.revision, request.revision);
    assert!(matches!(response.result,AnalysisResult::Search(hits) if hits[0].start==4));
    let mut stale = request;
    stale.revision = "wrong".into();
    assert!(
        process_request(&collagene_core::protocol::request::Request::encode(&stale).unwrap())
            .is_err()
    );
}

fn assembly_plan() -> (DesignPlan, Molecule) {
    let a = molecule("AAAACCCC", Topology::Linear);
    let b = molecule("CCCCGGGG", Topology::Linear);
    let target = molecule("AAAACCCCGGGG", Topology::Linear);
    let inputs = vec![
        PlanInput {
            id: "a".into(),
            revision: a.revision(),
            molecule: a,
        },
        PlanInput {
            id: "b".into(),
            revision: b.revision(),
            molecule: b,
        },
    ];
    (
        DesignPlan {
            standard: None,
            name: "independent assembly".into(),
            target_revision: target.revision(),
            inputs,
            steps: vec![Operation::Overlap {
                fragments: vec![Source::Input("a".into()), Source::Input("b".into())],
                overlaps: vec![4],
                circular: false,
            }],
        },
        target,
    )
}

#[test]
fn verifier_reconstructs_products_instead_of_believing_claims() {
    let (mut plan, target) = assembly_plan();
    let report = verify(&plan, &target);
    assert_eq!(report.status, Status::ConsistentUnderAssumptions);
    assert_eq!(
        report
            .product(0)
            .expect("selected step must reconstruct")
            .sequence,
        target.sequence
    );
    assert!(!report.exhaustive_alternative_search);
    if let Operation::Overlap { overlaps, .. } = &mut plan.steps[0] {
        overlaps[0] = 3;
    }
    assert_eq!(verify(&plan, &target).status, Status::Invalid);
}

#[test]
fn verifier_does_not_assume_overlapping_homology_arms_are_independent() {
    let a = molecule("AAAAAA", Topology::Linear);
    let target = molecule("AAAA", Topology::Circular);
    let plan = DesignPlan {
        standard: None,
        name: "overlapping arms".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "a".into(),
            revision: a.revision(),
            molecule: a,
        }],
        steps: vec![Operation::Overlap {
            fragments: vec![Source::Input("a".into()), Source::Input("a".into())],
            overlaps: vec![4, 4],
            circular: true,
        }],
    };
    assert_eq!(verify(&plan, &target).status, Status::Incomplete);
}

#[test]
fn verifier_rejects_stale_inputs_target_changes_and_cycles() {
    let (mut plan, target) = assembly_plan();
    plan.inputs[0].revision = "stale".into();
    assert_eq!(verify(&plan, &target).status, Status::Invalid);
    let (mut plan, target) = assembly_plan();
    plan.target_revision = "different".into();
    assert_eq!(verify(&plan, &target).status, Status::Invalid);
    let (mut plan, target) = assembly_plan();
    if let Operation::Overlap { fragments, .. } = &mut plan.steps[0] {
        fragments[0] = Source::Step(0);
    }
    assert_eq!(verify(&plan, &target).status, Status::Invalid);
}

#[test]
fn reports_escape_untrusted_names_and_bind_current_revisions() {
    let (mut plan, target) = assembly_plan();
    plan.name = "<script>alert(1)</script>".into();
    let html = export_protocol_html(&plan, &target);
    assert!(!html.contains("<script>"));
    assert!(html.contains("&lt;script&gt;"));
    assert!(html.contains(&target.revision()));
    assert!(html.contains(&collagene_core::planning::plan::DesignPlan::revision(&plan)));
}
