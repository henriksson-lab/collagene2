//! Independent literal top-call provenance fixtures; no target-guided feature
//! selection or physical activation is inferred by annotation correspondence.
use collagene_core::{
    Feature, Location, Molecule, Sequence, Topology,
    topo::{
        self as t,
        physical::{self as p, projection as v},
    },
    topo_annotations,
};

fn input(role: &str, dna: &str, topology: Topology, features: Vec<Feature>) -> t::Input {
    let mut molecule = Molecule::dna(role, dna, topology).unwrap();
    molecule.features = features;
    t::Input {
        role: role.into(),
        molecule,
    }
}
fn feature(label: &str, location: Location) -> Feature {
    Feature {
        kind: "misc_feature".into(),
        location,
        qualifiers: vec![("label".into(), Some(label.into()))],
        needs_review: false,
    }
}
fn strand(id: &str, dna: &str) -> p::Strand {
    p::Strand {
        role: id.into(),
        id: id.into(),
        sequence: Sequence::dna(dna).unwrap(),
        five_prime: p::FivePrime::Hydroxyl,
        three_prime: p::ThreePrime::Hydroxyl,
    }
}
fn correspondence(
    strand: usize,
    offset: usize,
    role: &str,
    start: usize,
    length: usize,
    reverse: bool,
) -> t::PreparedMap {
    t::PreparedMap {
        strand,
        strand_start: offset,
        length,
        origin: t::Origin::ReferenceCorrespondence(t::SourceInterval {
            source_role: role.into(),
            start,
            length,
            reverse,
        }),
    }
}
fn preparation(top: &str, bottom: &str, maps: Vec<t::PreparedMap>) -> t::PreparationEvidence {
    t::PreparationEvidence {
        substrates: p::Substrates {
            strands: vec![strand("top", top), strand("bottom", bottom)],
            pairing: vec![p::PairSpan {
                strand_a: 0,
                start_a: 0,
                strand_b: 1,
                start_b: 0,
                length: top.len(),
            }],
        },
        maps,
        not_represented: vec![],
        assumptions: vec!["Literal prepared correspondence; not biochemical proof".into()],
        retained_steps: vec![],
        layout: vec![],
        graph_origins: vec![],
    }
}
fn view(top: &str, bottom: &str) -> v::DuplexView {
    v::DuplexView {
        topology: v::ContourTopology::Linear,
        top: Sequence::dna(top).unwrap(),
        bottom: Sequence::dna(bottom).unwrap(),
        origin: p::Call {
            strand: 0,
            index: 0,
        },
        seed_reference_position: 0,
        source_maps: vec![v::PairedSourceMap {
            reference_start: 0,
            length: top.len(),
            top_strand: 0,
            top_start: 0,
            bottom_strand: 1,
            bottom_start: 0,
        }],
        nicks: vec![],
        ends: None,
        fully_known_terminal_chemistry: true,
        all_requested_steps_applied: true,
    }
}
fn project(
    inputs: &[t::Input],
    prep: &t::PreparationEvidence,
    view: &v::DuplexView,
) -> topo_annotations::AnnotationEvidence {
    let destination = Molecule::dna("destination", view.top.as_str(), Topology::Linear).unwrap();
    let before = serde_json::to_value((inputs, prep, view, &destination)).unwrap();
    let output = topo_annotations::project(inputs, prep, view, &destination).unwrap();
    assert_eq!(
        serde_json::to_value((inputs, prep, view, &destination)).unwrap(),
        before
    );
    output
}

#[test]
fn actual_origin_split_maps_project_only_top_parent_and_retain_unresolved_and_absent_features() {
    let inputs = vec![
        input(
            "template",
            "AAAACCCC",
            Topology::Circular,
            vec![
                feature(
                    "across origin",
                    Location::Join(vec![Location::range(6, 8), Location::range(0, 2)]),
                ),
                feature("outside", Location::range(2, 4)),
                feature("remote", Location::Unresolved("REMOTE.1:1..5".into())),
            ],
        ),
        input(
            "bottom_only",
            "TTGG",
            Topology::Linear,
            vec![feature("not a second parent", Location::range(0, 4))],
        ),
    ];
    let prep = preparation(
        "CCAA",
        "TTGG",
        vec![
            correspondence(0, 0, "template", 6, 4, false),
            correspondence(1, 0, "bottom_only", 0, 4, false),
        ],
    );
    let output = project(&inputs, &prep, &view("CCAA", "TTGG"));
    assert_eq!(output.features.len(), 1);
    assert_eq!(
        output.features[0]
            .location
            .extract(&Sequence::dna("CCAA").unwrap())
            .unwrap(),
        "CCAA"
    );
    assert_eq!(output.sources[0].maps.len(), 2);
    assert_eq!(
        (
            output.sources[0].maps[0].source_start,
            output.sources[0].maps[0].source_end,
            output.sources[0].maps[0].destination_start
        ),
        (0, 2, 2)
    );
    assert_eq!(
        (
            output.sources[0].maps[1].source_start,
            output.sources[0].maps[1].source_end,
            output.sources[0].maps[1].destination_start
        ),
        (6, 8, 0)
    );
    assert_eq!(output.sources[0].projected, 1);
    assert_eq!(output.sources[0].unmapped.len(), 2);
    assert!(
        output.sources[0]
            .unmapped
            .iter()
            .any(|u| u.source_index == 2 && u.reason.contains("Unresolved"))
    );
    assert!(output.sources[1].maps.is_empty());
    assert_eq!(output.sources[1].projected, 0);
    assert_eq!(
        output.sources[1].unmapped[0].feature.label(),
        "not a second parent"
    );
}

#[test]
fn reverse_correspondence_maps_actual_source_indices_without_copying_bottom_features() {
    let inputs = vec![
        input(
            "a",
            "ACGTTGCA",
            Topology::Linear,
            vec![feature("reverse", Location::range(2, 6))],
        ),
        input("b", "GTTG", Topology::Linear, vec![]),
    ];
    let prep = preparation(
        "CAAC",
        "GTTG",
        vec![
            correspondence(0, 0, "a", 2, 4, true),
            correspondence(1, 0, "b", 0, 4, false),
        ],
    );
    let output = project(&inputs, &prep, &view("CAAC", "GTTG"));
    assert_eq!(output.sources[0].maps[0].reverse, true);
    assert_eq!(output.features[0].location.spans(), vec![(0, 4, true)]);
    assert_eq!(
        output.features[0]
            .location
            .extract(&Sequence::dna("CAAC").unwrap())
            .unwrap(),
        "GTTG"
    );
}

#[test]
fn pcr_template_offsets_exclude_primer_and_new_tails_and_require_review() {
    for reverse in [false, true] {
        let (template, top, bottom, forward_tail, reverse_tail, tail) = if reverse {
            ("AGTC", "GGACTCC", "GGAGTCC", 2, 1, false)
        } else {
            ("ACGT", "GGACGTCCA", "TGGACGTCC", 2, 2, true)
        };
        let mut gene = feature("PCR template only", Location::range(0, 4));
        gene.kind = "CDS".into();
        gene.qualifiers.extend([
            ("translation".into(), Some("unverified translation".into())),
            ("codon_start".into(), Some("1".into())),
        ]);
        let inputs = vec![
            input("template", template, Topology::Linear, vec![gene]),
            input("bottom", bottom, Topology::Linear, vec![]),
        ];
        let copy_length = 4 + forward_tail + reverse_tail;
        let mut maps = vec![t::PreparedMap {
            strand: 0,
            strand_start: 0,
            length: copy_length,
            origin: t::Origin::ExactPcrCopy {
                material: 0,
                template_role: "template".into(),
                template_start: 0,
                template_length: 4,
                forward_tail_length: forward_tail,
                reverse_tail_length: reverse_tail,
                reverse,
            },
        }];
        if tail {
            maps.push(t::PreparedMap {
                strand: 0,
                strand_start: copy_length,
                length: 1,
                origin: t::Origin::SynthesizedTail {
                    material: 0,
                    base: collagene_core::tailing::TailingBase::A,
                },
            });
        }
        maps.push(correspondence(1, 0, "bottom", 0, bottom.len(), false));
        let output = project(&inputs, &preparation(top, bottom, maps), &view(top, bottom));
        let map = &output.sources[0].maps[0];
        assert_eq!(
            (
                map.source_start,
                map.source_end,
                map.destination_start,
                map.reverse
            ),
            (0, 4, if reverse { 1 } else { 2 }, reverse)
        );
        assert_eq!(output.features.len(), 1);
        assert!(output.features[0].needs_review);
        assert!(
            !output.features[0]
                .qualifiers
                .iter()
                .any(|(key, _)| key == "translation" || key == "codon_start")
        );
        assert_eq!(
            output.features[0]
                .location
                .extract(&Sequence::dna(top).unwrap())
                .unwrap(),
            template
        );
        assert!(
            output.sources[0]
                .warnings
                .iter()
                .any(|w| w.contains("primer tails"))
        );
    }
}

#[test]
fn repeated_source_copy_withholds_role_annotations_and_keeps_original_evidence() {
    let inputs = vec![
        input(
            "a",
            "AC",
            Topology::Linear,
            vec![feature("copied twice", Location::range(0, 2))],
        ),
        input("b", "GTGT", Topology::Linear, vec![]),
    ];
    let prep = preparation(
        "ACAC",
        "GTGT",
        vec![
            correspondence(0, 0, "a", 0, 2, false),
            correspondence(0, 2, "a", 0, 2, false),
            correspondence(1, 0, "b", 0, 4, false),
        ],
    );
    let output = project(&inputs, &prep, &view("ACAC", "GTGT"));
    assert!(output.features.is_empty());
    assert!(!output.sources[0].projection_supported);
    assert_eq!(output.sources[0].unmapped.len(), 1);
    assert_eq!(output.sources[0].unmapped[0].source_index, 0);
    assert_eq!(
        output.sources[0].unmapped[0].feature.qualifiers,
        inputs[0].molecule.features[0].qualifiers
    );
    assert!(output.sources[0].unmapped[0].reason.contains("repeated"));
}

#[test]
fn contradictory_destinations_call_ranges_and_resource_excess_fail_before_annotation_copies() {
    let inputs = vec![
        input(
            "a",
            "AC",
            Topology::Linear,
            vec![feature("a", Location::range(0, 2))],
        ),
        input("b", "GT", Topology::Linear, vec![]),
    ];
    let prep = preparation(
        "AC",
        "GT",
        vec![
            correspondence(0, 0, "a", 0, 2, false),
            correspondence(1, 0, "b", 0, 2, false),
        ],
    );
    let good = view("AC", "GT");
    let destination = Molecule::dna("bad", "TT", Topology::Linear).unwrap();
    assert!(topo_annotations::project(&inputs, &prep, &good, &destination).is_err());
    let destination = Molecule::dna("good", "AC", Topology::Linear).unwrap();
    let mut bad = good.clone();
    bad.source_maps[0].top_start = 1;
    assert!(topo_annotations::project(&inputs, &prep, &bad, &destination).is_err());
    let mut bad = prep.clone();
    bad.maps = vec![prep.maps[0].clone(); 4097];
    assert!(topo_annotations::project(&inputs, &bad, &good, &destination).is_err());
    let mut deep = Location::range(0, 1);
    for _ in 0..33 {
        deep = Location::Complement(Box::new(deep));
    }
    let mut deep_inputs = inputs.clone();
    deep_inputs[0].molecule.features[0].location = deep;
    assert!(topo_annotations::project(&deep_inputs, &prep, &good, &destination).is_err());
}
