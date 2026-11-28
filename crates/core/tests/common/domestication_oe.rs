use collagene_core::{
    cloning::{
        golden_gate::domestication::{self as design, preparation as physical},
        overlap::extension::primers as overlap_extension_primer_design,
        restriction::catalog as catalog_digest,
    },
    planning::plan::PlanInput,
    sequence::codon_design,
};
use collagene_model::{Molecule, Topology};
pub const ORIGINAL: &str = "AGTCCGATGCTACGTAGAATTCCAGACCTGTAGATCGTCA";
pub const INTENDED: &str = "AGTCCGATGCTACGTAGACTTCCAGACCTGTAGATCGTCA";
pub fn fixture(topology: Topology) -> (PlanInput, physical::Request) {
    let source = Molecule::dna("Original40 EcoRI source", ORIGINAL, topology).unwrap();
    let design = design::Request {
        source_revision: source.revision(),
        selection: catalog_digest::select(
            &catalog_digest::builtin_catalog().unwrap(),
            &["ecori".into()],
        )
        .unwrap(),
        scope: vec![design::SourceRange { start: 0, end: 40 }],
        retained_sites: vec![],
        coding_features: vec![],
        protected: vec![],
        options: codon_design::DesignOptions::default(),
        limits: design::Limits::default(),
    };
    let patches = vec![design::Patch {
        position: 18,
        expected: "A".into(),
        replacement: "C".into(),
    }];
    let selection = design::Selection {
        context_revision: design::propose_revision(&source, &design).unwrap(),
        candidate_id: design::candidate_id(&design, &patches).unwrap(),
        patches,
    };
    let options = overlap_extension_primer_design::Options {
        primer_options: collagene_core::cloning::pcr::primers::DesignOptions {
            min_length: 8,
            max_length: 8,
            min_gc_fraction: 0.,
            max_gc_fraction: 1.,
            max_tm_difference_c: 100.,
            max_homopolymer: 200,
            max_three_prime_complementarity: 200,
            ..Default::default()
        },
        homology_flank: 8,
        min_paired_bases: 8,
        ..Default::default()
    };
    (
        PlanInput {
            id: physical::INPUT_ID.into(),
            revision: source.revision(),
            molecule: source,
        },
        physical::Request {
            origin: None,
            design,
            selection,
            route: physical::Route::OverlapExtension { options },
        },
    )
}
