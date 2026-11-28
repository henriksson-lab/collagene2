//! NCBI gc.prt v4.6 table data, in TCAG codon order. See docs/GENETIC_CODE_CATALOG.md.
//! The Starts row marks both initiation (M) and possible termination (*).
use super::CodeDescriptor;
use collagene_model::Error;

#[derive(Clone, Copy)]
pub(super) struct Code {
    pub(super) id: u8,
    pub(super) name: &'static str,
    pub(super) amino_acids: &'static [u8; 64],
    starts: &'static [u8; 64],
}
impl Code {
    pub(super) fn lookup(id: u8) -> Result<Self, Error> {
        CODES
            .iter()
            .find(|code| code.id == id)
            .copied()
            .ok_or_else(|| {
                Error::Unsupported(format!(
                    "genetic code {id} is not in the NCBI v4.6 catalogue"
                ))
            })
    }
    pub(super) fn is_start(self, index: usize) -> bool {
        self.starts[index] == b'M'
    }
    pub(super) fn context_stop(self, index: usize) -> bool {
        self.starts[index] == b'*' && self.amino_acids[index] != b'*'
    }
    pub(super) fn has_context_stops(self) -> bool {
        (0..64).any(|index| self.context_stop(index))
    }
    pub(super) fn descriptor(self) -> CodeDescriptor {
        CodeDescriptor {
            id: self.id,
            name: self.name.into(),
            start_codons: (0..64)
                .filter(|&i| self.is_start(i))
                .map(index_codon)
                .collect(),
            context_dependent_stop_codons: (0..64)
                .filter(|&i| self.context_stop(i))
                .map(index_codon)
                .collect(),
        }
    }
}
fn index_codon(index: usize) -> String {
    let bases = b"TCAG";
    [bases[index / 16], bases[index / 4 % 4], bases[index % 4]]
        .into_iter()
        .map(char::from)
        .collect()
}
pub(super) fn descriptors() -> Vec<CodeDescriptor> {
    CODES.iter().map(|code| code.descriptor()).collect()
}
const CODES: &[Code] = &[
    Code {
        id: 1,
        name: "Standard",
        amino_acids: b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"---M------**--*----M---------------M----------------------------",
    },
    Code {
        id: 2,
        name: "Vertebrate mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
        starts: b"----------**--------------------MMMM----------**---M------------",
    },
    Code {
        id: 3,
        name: "Yeast Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**----------------------MM---------------M------------",
    },
    Code {
        id: 4,
        name: "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--MM------**-------M------------MMMM---------------M------------",
    },
    Code {
        id: 5,
        name: "Invertebrate Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
        starts: b"---M------**--------------------MMMM---------------M------------",
    },
    Code {
        id: 6,
        name: "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear",
        amino_acids: b"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--------------*--------------------M----------------------------",
    },
    Code {
        id: 9,
        name: "Echinoderm Mitochondrial; Flatworm Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        starts: b"----------**-----------------------M---------------M------------",
    },
    Code {
        id: 10,
        name: "Euplotid Nuclear",
        amino_acids: b"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**-----------------------M----------------------------",
    },
    Code {
        id: 11,
        name: "Bacterial, archaeal and plant plastid",
        amino_acids: b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"---M------**--*----M------------MMMM---------------M------------",
    },
    Code {
        id: 12,
        name: "Alternative Yeast Nuclear",
        amino_acids: b"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**--*----M---------------M----------------------------",
    },
    Code {
        id: 13,
        name: "Ascidian Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
        starts: b"---M------**----------------------MM---------------M------------",
    },
    Code {
        id: 14,
        name: "Alternative Flatworm Mitochondrial",
        amino_acids: b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        starts: b"-----------*-----------------------M----------------------------",
    },
    Code {
        id: 15,
        name: "Blepharisma Macronuclear",
        amino_acids: b"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------*---*--------------------M----------------------------",
    },
    Code {
        id: 16,
        name: "Chlorophycean Mitochondrial",
        amino_acids: b"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------*---*--------------------M----------------------------",
    },
    Code {
        id: 21,
        name: "Trematode Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        starts: b"----------**-----------------------M---------------M------------",
    },
    Code {
        id: 22,
        name: "Scenedesmus obliquus Mitochondrial",
        amino_acids: b"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"------*---*---*--------------------M----------------------------",
    },
    Code {
        id: 23,
        name: "Thraustochytrium Mitochondrial",
        amino_acids: b"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--*-------**--*-----------------M--M---------------M------------",
    },
    Code {
        id: 24,
        name: "Rhabdopleuridae Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        starts: b"---M------**-------M---------------M---------------M------------",
    },
    Code {
        id: 25,
        name: "Candidate Division SR1 and Gracilibacteria",
        amino_acids: b"FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"---M------**-----------------------M---------------M------------",
    },
    Code {
        id: 26,
        name: "Pachysolen tannophilus Nuclear",
        amino_acids: b"FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**--*----M---------------M----------------------------",
    },
    Code {
        id: 27,
        name: "Karyorelict Nuclear",
        amino_acids: b"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--------------*--------------------M----------------------------",
    },
    Code {
        id: 28,
        name: "Condylostoma Nuclear",
        amino_acids: b"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**--*--------------------M----------------------------",
    },
    Code {
        id: 29,
        name: "Mesodinium Nuclear",
        amino_acids: b"FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--------------*--------------------M----------------------------",
    },
    Code {
        id: 30,
        name: "Peritrich Nuclear",
        amino_acids: b"FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--------------*--------------------M----------------------------",
    },
    Code {
        id: 31,
        name: "Blastocrithidia Nuclear",
        amino_acids: b"FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**-----------------------M----------------------------",
    },
    Code {
        id: 32,
        name: "Balanophoraceae Plastid",
        amino_acids: b"FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"---M------*---*----M------------MMMM---------------M------------",
    },
    Code {
        id: 33,
        name: "Cephalodiscidae Mitochondrial",
        amino_acids: b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        starts: b"---M-------*-------M---------------M---------------M------------",
    },
];
