//! Version-pinned CC0 human V-gene prefixes with provider-supplied IMGT labels.
//! These are germline reference facts, not complete rearranged domains, trained
//! classification profiles, species predictions or automatic numbering results.
use super::*;
use collagene_model::PROTEIN_RESIDUES;
use std::collections::BTreeSet;
pub const REFERENCE_ENGINE: &str = "antibody-reference-corpus-1";
const DATA: &str = include_str!("../data/ogrdb-human-v/references.json");
data! {pub struct Excluded {pub allele_ref:String,pub reason:String,pub sequence_type:String}}
data! {pub struct Source {
    pub set_id:String,pub set_ref:String,pub release_version:u32,pub release_date:String,
    pub species:String,pub taxonomy_id:String,pub locus:String,pub author:String,
    pub url:String,pub sha256:String,pub original_records:usize,pub retained_records:usize,
    pub excluded:Vec<Excluded>,
}}
data! {pub struct Cdr {
    pub kind:u8,
    /// Zero-based half-open offsets in the ungapped V-prefix protein.
    pub start:usize,pub end:usize,
}}
data! {pub struct Reference {
    pub source:usize,pub allele_ref:String,pub allele_id:String,pub label:String,
    pub release_version:u32,
    /// Exactly 104 codons. A whole gap codon is represented by three dots.
    pub aligned_nucleotides:String,
    /// Exactly 104 columns; a '-' denotes a missing position, never a residue.
    pub aligned_protein:String,
    /// Explicit provider labels, not labels inferred from a PDB author's IDs.
    pub alignment_labels:Vec<String>,pub cdrs:Vec<Cdr>,
    /// Retained atypical functional annotations. Do not silently repair these
    /// residues or accept this reference as a canonical conserved-anchor model.
    pub non_cysteine_anchors:Vec<u16>,
}}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Corpus {
    version: u16,
    license: String,
    license_url: String,
    scope: String,
    sources: Vec<Source>,
    entries: Vec<Reference>,
}
impl Corpus {
    pub fn sources(&self) -> &[Source] {
        &self.sources
    }
    pub fn entries(&self) -> &[Reference] {
        &self.entries
    }
    pub fn revision(&self) -> Result<String> {
        self.validate()?;
        hash(
            b"antibody-reference-corpus-1\0",
            &(REFERENCE_ENGINE, self),
            1 << 20,
        )
    }
    fn validate(&self) -> Result<()> {
        if self.version != 1
            || self.license != "CC0-1.0"
            || self.license_url
                != "https://ogrdb.airr-community.org/render_page/licensing_statement.html"
            || self.scope != "human-functional-V-IMGT-1-through-104"
            || self.sources.len() != 3
            || self.entries.len() != 354
        {
            return Err(Error::Invalid("unsupported embedded reference corpus"));
        }
        let expected = [
            ("OGRDB:G00060", 10, "IGH", 194),
            ("OGRDB:G00061", 5, "IGK", 74),
            ("OGRDB:G00062", 4, "IGL", 86),
        ];
        let mut identities = BTreeSet::new();
        let mut counts = [0; 3];
        for (source, &(id, version, locus, retained)) in self.sources.iter().zip(&expected) {
            if source.set_id != id
                || source.release_version != version
                || source.locus != locus
                || source.retained_records != retained
                || source.species != "Homo sapiens"
                || source.taxonomy_id != "NCBITAXON:9606"
                || source.author.trim().is_empty()
                || source.sha256.len() != 64
                || !source
                    .sha256
                    .bytes()
                    .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
                || source.retained_records + source.excluded.len() != source.original_records
            {
                return Err(Error::Invalid("reference source metadata"));
            }
            for excluded in &source.excluded {
                if excluded.reason != "not-functional-V"
                    || !identities.insert(excluded.allele_ref.as_str())
                {
                    return Err(Error::Invalid("reference exclusion ledger"));
                }
            }
        }
        for entry in &self.entries {
            if entry.source >= 3
                || entry.release_version == 0
                || entry.label.is_empty()
                || entry.allele_id.is_empty()
                || !identities.insert(entry.allele_ref.as_str())
                || entry.aligned_nucleotides.len() != 312
                || entry.aligned_protein.len() != 104
                || entry.alignment_labels != (1..=104).map(|n| n.to_string()).collect::<Vec<_>>()
            {
                return Err(Error::Invalid("reference identity/columns"));
            }
            counts[entry.source] += 1;
            let mut labels = vec![];
            for (column, (&aa, dna)) in entry
                .aligned_protein
                .as_bytes()
                .iter()
                .zip(entry.aligned_nucleotides.as_bytes().chunks_exact(3))
                .enumerate()
            {
                if aa == b'-' {
                    if dna != b"..." {
                        return Err(Error::Invalid("reference gap codon"));
                    }
                } else {
                    if !PROTEIN_RESIDUES.contains(&aa) || !dna.iter().all(|b| b"ACGT".contains(b)) {
                        return Err(Error::Invalid("reference alphabet/codon"));
                    }
                    labels.push(column + 1);
                }
            }
            if entry.non_cysteine_anchors
                != [23_u16, 104]
                    .into_iter()
                    .filter(|&n| entry.aligned_protein.as_bytes()[usize::from(n) - 1] != b'C')
                    .collect::<Vec<_>>()
            {
                return Err(Error::Invalid("reference atypical anchor ledger"));
            }
            if entry.cdrs.len() != 2 {
                return Err(Error::Invalid("reference CDR delineations"));
            }
            for (cdr, &(kind, low, high)) in entry.cdrs.iter().zip(&[(1, 27, 38), (2, 56, 65)]) {
                let selected = labels
                    .iter()
                    .enumerate()
                    .filter_map(|(i, &n)| (low..=high).contains(&n).then_some(i))
                    .collect::<Vec<_>>();
                if cdr.kind != kind
                    || selected.first() != Some(&cdr.start)
                    || selected.last().map(|n| n + 1) != Some(cdr.end)
                {
                    return Err(Error::Invalid("reference CDR/label disagreement"));
                }
            }
        }
        if counts != [194, 74, 86] {
            return Err(Error::Invalid("reference source counts"));
        }
        Ok(())
    }
}
/// Load the finite embedded corpus. No network request, executable or database
/// installation is performed. Nucleotide/protein translations are independently
/// checked against the existing scientific core in the corpus integration test.
pub fn human_v_references() -> Result<Corpus> {
    if DATA.len() > 1 << 20 {
        return Err(Error::Limit("embedded reference bytes"));
    }
    let corpus: Corpus =
        serde_json::from_str(DATA).map_err(|_| Error::Invalid("embedded reference JSON"))?;
    corpus.validate()?;
    Ok(corpus)
}
