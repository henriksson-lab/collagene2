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
pub mod transfer;
pub mod verification;

pub use model::{Alphabet, Error, Feature, Location, Molecule, Sequence, Topology};

pub mod embl;
pub mod snapgene;
pub mod traces;

pub mod svg;

pub mod annotations;

pub mod translation;

pub mod alignment;

pub mod editing;

pub mod primers;

pub mod pdf;

pub mod gel;

pub mod mutagenesis;

pub mod golden_gate;

pub mod tracks;

pub mod workspace;

pub mod history;

pub mod oligo_formats;

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

pub mod restriction_history;

/// Explicit source-backed physical TOPO preparation and reconstruction.
pub mod topo;

pub mod read_alignment_reference;
pub mod read_alignment_store;

pub mod read_alignment_archive;

pub mod read_alignment_session;
