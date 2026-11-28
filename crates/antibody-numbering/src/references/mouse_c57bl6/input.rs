//! The pinned external OGRDB representation, deliberately not an internal format.
//! Only fields needed for this curation are decoded; untouched original bytes
//! and their exact identities remain separate acquisition evidence.
use crate::{Error, Result};
use serde::Deserialize;

pub(super) const MAX_JSON: usize = 512 * 1024;
pub(super) const MAX_FASTA: usize = 64 * 1024;
#[derive(Clone, Deserialize)]
pub(super) struct Envelope {
    #[serde(rename = "GermlineSet")]
    pub sets: Vec<RawSet>,
}
#[derive(Clone, Deserialize)]
pub(super) struct Taxon {
    pub id: String,
    pub label: String,
}
#[derive(Clone, Deserialize)]
pub(super) struct RawSet {
    pub germline_set_id: String,
    pub germline_set_ref: String,
    pub release_version: serde_json::Number,
    pub release_date: String,
    pub species: Taxon,
    pub species_subgroup: String,
    pub species_subgroup_type: String,
    pub locus: String,
    pub author: String,
    pub allele_descriptions: Vec<Allele>,
}
#[derive(Clone, Deserialize)]
pub(super) struct Allele {
    pub allele_description_id: String,
    pub allele_description_ref: String,
    pub release_version: serde_json::Number,
    pub label: String,
    pub sequence: String,
    pub coding_sequence: String,
    pub locus: String,
    pub species: Taxon,
    pub species_subgroup: String,
    pub sequence_type: String,
    pub functional: bool,
    pub gene_start: Option<usize>,
    pub gene_end: Option<usize>,
    pub j_codon_frame: Option<u8>,
    pub j_cdr3_end: Option<usize>,
    pub v_gene_delineations: Option<Vec<Delineation>>,
}
#[derive(Clone, Deserialize)]
pub(super) struct Delineation {
    pub delineation_scheme: String,
    pub aligned_sequence: Option<String>,
    pub alignment_labels: Option<Vec<String>>,
    pub cdr1_start: usize,
    pub cdr1_end: usize,
    pub cdr2_start: usize,
    pub cdr2_end: usize,
    pub fwr3_end: usize,
}
#[derive(Clone)]
pub(super) struct GappedRecord {
    pub label: String,
    pub sequence: String,
}

pub(super) fn version(value: &serde_json::Number) -> Result<u32> {
    let n = value
        .as_f64()
        .filter(|n| n.is_finite() && *n >= 1. && *n <= 1_000. && n.fract() == 0.)
        .ok_or(Error::Invalid("mouse source release version"))?;
    Ok(n as u32)
}
pub(super) fn decode(json: &[u8], fasta: &[u8]) -> Result<(RawSet, Vec<GappedRecord>)> {
    if json.len() > MAX_JSON || fasta.len() > MAX_FASTA {
        return Err(Error::Limit("mouse source bytes"));
    }
    let mut envelope: Envelope =
        serde_json::from_slice(json).map_err(|_| Error::Invalid("mouse source JSON"))?;
    if envelope.sets.len() != 1 {
        return Err(Error::Invalid("one mouse source set required"));
    }
    let text = std::str::from_utf8(fasta).map_err(|_| Error::Invalid("mouse FASTA UTF-8"))?;
    let mut records = Vec::<GappedRecord>::new();
    for line in text.lines() {
        if let Some(label) = line.strip_prefix('>') {
            if records.len() >= 114
                || label.is_empty()
                || label.len() > 128
                || !label.bytes().all(|b| b.is_ascii_graphic())
            {
                return Err(Error::Invalid("mouse FASTA label/count"));
            }
            records.push(GappedRecord {
                label: label.into(),
                sequence: String::new(),
            });
        } else {
            let current = records
                .last_mut()
                .ok_or(Error::Invalid("mouse FASTA sequence before header"))?;
            if line.is_empty()
                || !line.bytes().all(|b| b"ACGT.".contains(&b))
                || line.len() > 2_048_usize.saturating_sub(current.sequence.len())
            {
                return Err(Error::Invalid("mouse FASTA alphabet/length"));
            }
            current.sequence.push_str(line);
        }
    }
    if records.len() != 114 || records.iter().any(|r| r.sequence.is_empty()) {
        return Err(Error::Invalid("mouse FASTA record inventory"));
    }
    Ok((envelope.sets.remove(0), records))
}
