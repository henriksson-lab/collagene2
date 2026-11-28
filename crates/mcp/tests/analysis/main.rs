//! MCP analysis integration tests, grouped by biological responsibility.

mod alignment_tools;
mod antibody_tools;
mod base_editing_tools;
mod codon_tools;
mod contig_tools;
mod crispr_tools;
mod design_checks;
mod dotplot_tools;
mod genome_tools;
mod ivt_linearization_tools;
mod ivt_tools;
mod methylation_tools;
mod protein_liability_tools;
mod read_alignment_tools;
mod rna_chemistry_tools;
mod scientific_tools;
mod sequence_similarity_tools;
mod structure_tools;
mod variant_tools;
#[path = "vector_class/mod.rs"]
mod vector_class_tools;

#[path = "causal_qc/mod.rs"]
mod causal_qc_tools;
