//! Independent literal coding oracles for supplied-consensus comparisons.
//! These are input observations, never native edit/assembly ancestry receipts.
use collagene_model::{Feature, Location, Molecule, Topology};

pub(super) struct CodingLiteral {
    pub name: &'static str,
    pub design: &'static str,
    pub consensus: &'static str,
    pub location: Location,
    pub topology: Topology,
    /// Complete biological coding strings expected under the selected path.
    pub design_coding: &'static str,
    pub consensus_coding: &'static str,
    /// Internal codon ordinal; None means no new definite internal stop.
    pub new_stop: Option<usize>,
}
impl CodingLiteral {
    pub fn molecules(&self) -> (Molecule, Molecule) {
        let mut design = Molecule::dna(self.name, self.design, self.topology).unwrap();
        design.features.push(Feature {
            kind: "CDS".into(),
            location: self.location.clone(),
            qualifiers: vec![("codon_start".into(), Some("1".into()))],
            needs_review: false,
        });
        let consensus = Molecule::dna(
            "Supplied unannotated consensus",
            self.consensus,
            self.topology,
        )
        .unwrap();
        (design, consensus)
    }
}

pub(super) fn coding_literals() -> Vec<CodingLiteral> {
    vec![
        CodingLiteral {
            name: "Synonymous GAA to GAG at codon1",
            design: "ATGGAACTGTAA",
            consensus: "ATGGAGCTGTAA",
            location: Location::range(0, 12),
            topology: Topology::Linear,
            design_coding: "ATGGAACTGTAA",
            consensus_coding: "ATGGAGCTGTAA",
            new_stop: None,
        },
        CodingLiteral {
            name: "CAA to definite internal TAA at codon1",
            design: "ATGCAACCCTAA",
            consensus: "ATGTAACCCTAA",
            location: Location::range(0, 12),
            topology: Topology::Linear,
            design_coding: "ATGCAACCCTAA",
            consensus_coding: "ATGTAACCCTAA",
            new_stop: Some(1),
        },
        CodingLiteral {
            name: "One inserted A produces distant new TAA at codon4",
            design: "ATGAAACCAGGTAACGAATAA",
            consensus: "ATGAAAACCAGGTAACGAATAA",
            location: Location::range(0, 21),
            topology: Topology::Linear,
            design_coding: "ATGAAACCAGGTAACGAATAA",
            consensus_coding: "ATGAAAACCAGGTAACGAATAA",
            new_stop: Some(4),
        },
        CodingLiteral {
            name: "Three inserted GCT calls preserve downstream register",
            design: "ATGCAACCCTAA",
            consensus: "ATGCAAGCTCCCTAA",
            location: Location::range(0, 12),
            topology: Topology::Linear,
            design_coding: "ATGCAACCCTAA",
            consensus_coding: "ATGCAAGCTCCCTAA",
            new_stop: None,
        },
        CodingLiteral {
            name: "Reverse CDS literal CAA to TAA",
            design: "TTAGGGTTGCAT",
            consensus: "TTAGGGTTACAT",
            location: Location::Complement(Box::new(Location::range(0, 12))),
            topology: Topology::Linear,
            design_coding: "ATGCAACCCTAA",
            consensus_coding: "ATGTAACCCTAA",
            new_stop: Some(1),
        },
        CodingLiteral {
            name: "Joined CDS excludes the literal GG intron",
            design: "ATGCAAGGCCCTAA",
            consensus: "ATGTAAGGCCCTAA",
            location: Location::Join(vec![Location::range(0, 6), Location::range(8, 14)]),
            topology: Topology::Linear,
            design_coding: "ATGCAACCCTAA",
            consensus_coding: "ATGTAACCCTAA",
            new_stop: Some(1),
        },
        CodingLiteral {
            name: "Intron insertion alone does not shift coding register",
            design: "ATGCAAGGCCCTAA",
            consensus: "ATGCAAGAGCCCTAA",
            location: Location::Join(vec![Location::range(0, 6), Location::range(8, 14)]),
            topology: Topology::Linear,
            design_coding: "ATGCAACCCTAA",
            consensus_coding: "ATGCAACCCTAA",
            new_stop: None,
        },
        CodingLiteral {
            name: "Origin-spanning CDS uses 9..15 then0..6",
            design: "CCCTAACGTATGCAA",
            consensus: "CCCTAACGTATGTAA",
            location: Location::Join(vec![Location::range(9, 15), Location::range(0, 6)]),
            topology: Topology::Circular,
            design_coding: "ATGCAACCCTAA",
            consensus_coding: "ATGTAACCCTAA",
            new_stop: Some(1),
        },
        CodingLiteral {
            name: "TAN has both tyrosine and stop expansions",
            design: "ATGCAACCCTAA",
            consensus: "ATGTANCCCTAA",
            location: Location::range(0, 12),
            topology: Topology::Linear,
            design_coding: "ATGCAACCCTAA",
            consensus_coding: "ATGTANCCCTAA",
            new_stop: None,
        },
        CodingLiteral {
            name: "TAR has only stop expansions under table1",
            design: "ATGCAACCCTAA",
            consensus: "ATGTARCCCTAA",
            location: Location::range(0, 12),
            topology: Topology::Linear,
            design_coding: "ATGCAACCCTAA",
            consensus_coding: "ATGTARCCCTAA",
            new_stop: Some(1),
        },
        CodingLiteral {
            name: "Compensating insertion and deletion restore downstream register",
            design: "ATGAAACCAGGTAACGAATAA",
            consensus: "ATGAAAACCAGTAACGAATAA",
            location: Location::range(0, 21),
            topology: Topology::Linear,
            design_coding: "ATGAAACCAGGTAACGAATAA",
            consensus_coding: "ATGAAAACCAGTAACGAATAA",
            new_stop: None,
        },
        CodingLiteral {
            name: "Reverse joined CDS excludes CC and preserves biological order",
            design: "TTAGGGCCTTGCAT",
            consensus: "TTAGGGCCTTACAT",
            location: Location::Complement(Box::new(Location::Join(vec![
                Location::range(0, 6),
                Location::range(8, 14),
            ]))),
            topology: Topology::Linear,
            design_coding: "ATGCAACCCTAA",
            consensus_coding: "ATGTAACCCTAA",
            new_stop: Some(1),
        },
        CodingLiteral {
            name: "A query insertion at stored origin is internal to the circular CDS",
            design: "CCCTAACGTATGCAA",
            consensus: "ACCCTAACGTATGCAA",
            location: Location::Join(vec![Location::range(9, 15), Location::range(0, 6)]),
            topology: Topology::Circular,
            design_coding: "ATGCAACCCTAA",
            consensus_coding: "ATGCAAACCCTAA",
            new_stop: None,
        },
    ]
}

/// Literal table1 translations, including terminal stops. No translation or
/// alignment algorithm is used to manufacture the expected test result.
pub(super) fn expected_protein(coding: &str) -> &'static str {
    match coding {
        "ATGGAACTGTAA" | "ATGGAGCTGTAA" => "MEL*",
        "ATGCAACCCTAA" => "MQP*",
        "ATGTAACCCTAA" | "ATGTARCCCTAA" => "M*P*",
        "ATGTANCCCTAA" => "MXP*",
        "ATGCAAGCTCCCTAA" => "MQAP*",
        "ATGCAAACCCTAA" => "MQTL",
        "ATGAAACCAGGTAACGAATAA" => "MKPGNE*",
        "ATGAAAACCAGGTAACGAATAA" => "MKTR*RI",
        "ATGAAAACCAGTAACGAATAA" => "MKTSNE*",
        _ => panic!("missing independently authored coding oracle: {coding}"),
    }
}
