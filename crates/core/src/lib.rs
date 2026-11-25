//! Portable scientific model and independently implemented cloning algorithms.
//! No UI, filesystem, subprocess, network or async runtime dependencies.

pub mod analysis;
pub mod catalog_digest;
pub mod chemistry;
pub mod formats;
pub mod model;
pub mod protocol;
pub mod recombination;
pub mod restriction_inspection;
pub mod restriction_selection;
pub mod transfer;
pub mod verification;

pub use model::{Alphabet, Error, Feature, Location, Molecule, Sequence, Topology};

pub mod embl;
pub mod snapgene;
pub mod traces;

pub mod svg;

pub mod annotation_batch;
pub mod annotations;

pub mod translation;
pub mod translation_reference;

pub mod alignment;
pub mod pairwise_review;

pub mod editing;

pub mod primers;

pub mod pdf;

pub mod gel;

pub mod mutagenesis;

pub mod golden_gate;

pub mod tracks;

pub mod workspace;

pub mod history;

pub mod mismatch_pcr_outputs;
pub mod oligo_formats;
pub mod oligo_orders;
pub mod pcr_outputs;

pub mod trace_validation;

pub mod feature_library;

pub mod annotation_formats;

pub mod legacy_formats;

pub mod search;

pub mod alignment_formats;

pub mod view_tracks;

pub mod protein_properties;

pub mod tailing;

pub mod overlap_profiles;

pub mod run_lineage;

pub mod codon_design;

pub mod dotplot;

pub mod contig;

pub mod contig_history;

pub mod methylation;

pub mod modified_duplex;

pub mod edit_template;

pub mod crispr;

pub mod genome;

pub mod genome_crispr;

pub mod genome_session;

pub mod digest_setups;
pub mod restriction_history;

/// Explicit source-backed physical TOPO preparation and reconstruction.
pub mod topo;
pub mod topo_annotations;

pub mod indexed_bam;
pub mod read_alignment_reference;
pub mod read_alignment_store;

pub mod read_alignment_archive;

pub mod read_alignment_session;

pub mod topo_worker;

pub mod variant_reference;

pub mod variants_archive;

pub mod variants_session;

pub mod assembly_recipe;

pub mod primer_reuse;

pub mod ligation_recipe;

pub mod assembly_junctions;

pub mod overlap_primer_design;

pub mod assembly_intent;
pub mod combinatorial_batch;
pub mod combinatorial_design;
pub mod combinatorial_export;
pub mod combinatorial_release;

pub mod digest_cut_selection;

pub mod cloning_outputs;
pub mod plan_oligos;

pub mod golden_gate_primer_design;

pub mod mutation_intent;

pub mod read_checks;

pub mod read_review;

pub mod assembly_junction_edit;
pub mod assembly_library_history;

pub mod selective_junction_edit;

pub mod assembly_diagnostics;

pub mod domestication;
pub mod domestication_preparation;
pub mod golden_gate_fusion_design;
pub mod golden_gate_fusion_rank;
pub mod golden_gate_junction_edit;
pub mod ligation_fidelity;
pub mod ligation_profiles;
pub mod topo_preparation;

pub mod prepared_golden_gate;

pub mod golden_gate_reaction;

pub mod pairwise_context;
pub mod pairwise_export;
pub mod pairwise_history;
pub mod pairwise_search;

pub mod reference_transfer;

pub mod enzyme_sets;

pub mod topo_nick_sealing;

pub mod topo_processed_product;

mod protocol_legacy_v64;

pub mod processed_run_evidence;

pub mod topo_actions;

pub mod physical_run_evidence;

/// Dependency-preserving manual edits to declarative cloning plans.
pub mod plan_editing;

/// Source-bound multi-parent reference recombination and annotation evidence.
pub mod multisite_recombination;

pub mod reference_run_evidence;
