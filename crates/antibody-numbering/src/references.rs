//! Version-pinned CC0 species/strain-scoped V-prefix and J/FR4 reference facts.
//! These are germline reference facts, not complete rearranged domains, trained
//! classification profiles, species predictions or automatic numbering results.
use super::*;
use collagene_model::PROTEIN_RESIDUES;
use std::collections::BTreeSet;
pub const REFERENCE_ENGINE: &str = "antibody-reference-corpus-3";
#[cfg(any(test, feature = "reference-curation"))]
pub mod mouse_c57bl6;
data! {
    pub enum Species {
        Human,
        MouseC57Bl6
    }
}
data! {
    pub struct AtypicalAnchor {
        pub allele_ref: String,
        pub imgt_column: u16,
        pub codon: String,
        pub residue: char,
        pub canonical_residue: char
    }
}
data! {
    pub enum SourcePolicy {
        HumanProviderLabeled,
        MousePinnedColumns {
            gapped_sha256: String,
            v_label_evidence: String,
            #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, AtypicalAnchor, 4>")]
            atypical_anchors: Vec<AtypicalAnchor>,
            j_phase_policy: String,
        },
    }
}
const HUMAN_DATA: &str = include_str!("../data/ogrdb-human-vj/references.json");
const MOUSE_DATA: &str = include_str!("../data/ogrdb-mouse-c57bl6-igh-v5/references.json");
data! {
    pub struct Excluded {
        pub allele_ref: String,
        pub reason: String,
        pub sequence_type: String
    }
}
data! {
    pub struct Source {
        pub set_id: String,
        pub set_ref: String,
        pub release_version: u32,
        pub release_date: String,
        pub species: String,
        pub taxonomy_id: String,
        pub locus: String,
        pub author: String,
        pub url: String,
        pub sha256: String,
        pub original_records: usize,
        pub v_records: usize,
        pub j_records: usize,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Excluded, 64>")]
        pub excluded: Vec<Excluded>,
    }
}
data! {
    pub struct Cdr {
        pub kind: u8,
        /// Zero-based half-open offsets in the ungapped V-prefix protein.
        pub start: usize,
        pub end: usize,
    }
}
data! {
    pub struct VReference {
        pub source: usize,
        pub allele_ref: String,
        pub allele_id: String,
        pub label: String,
        pub release_version: u32,
        /// Exactly 104 codons. A whole gap codon is represented by three dots.
        pub aligned_nucleotides: String,
        /// Exactly 104 columns; a '-' denotes a missing position, never a residue.
        pub aligned_protein: String,
        /// Column labels with the origin declared by Corpus::source_policy.
        /// Mouse labels are authored IMGT-column ordinals, not JSON provider labels.
        pub alignment_labels: Vec<String>,
        pub cdrs: Vec<Cdr>,
        /// Retained atypical functional annotations. Do not silently repair these
        /// residues or accept this reference as a canonical conserved-anchor model.
        pub non_cysteine_anchors: Vec<u16>,
    }
}
data! {
    pub struct JReference {
        pub source: usize,
        pub allele_ref: String,
        pub allele_id: String,
        pub label: String,
        pub release_version: u32,
        pub original_sequence: String,
        pub coding_nucleotides: String,
        /// Original provider gene interval: one-based inclusive, in original_sequence.
        pub gene_start: usize,
        pub gene_end: usize,
        /// Raw OGRDB values, retained without claiming generic AIRR interpretation.
        pub provider_codon_frame: u8,
        pub provider_j_cdr3_end: usize,
        /// Zero-based coding-sequence offset verified for these pinned releases.
        pub anchor_offset: usize,
        pub fr4_nucleotides: String,
        pub fr4_protein: String,
        /// Authored IMGT-118 sequential rule, not provider-supplied position labels.
        pub fr4_positions: Vec<String>,
        /// Incomplete splice-junction codon; never translated or numbered as FR4.
        pub trailing_nucleotides: String,
    }
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Corpus {
    version: u16,
    species: Species,
    source_policy: SourcePolicy,
    license: String,
    license_url: String,
    scope: String,
    j_label_rule: String,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Source, 3>")]
    sources: Vec<Source>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, VReference, 354>")]
    v_entries: Vec<VReference>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, JReference, 23>")]
    j_entries: Vec<JReference>,
}
impl Corpus {
    pub fn species(&self) -> &Species {
        &self.species
    }
    pub fn source_policy(&self) -> &SourcePolicy {
        &self.source_policy
    }
    /// Authored conditional FR4 rule, not provider-supplied J position labels.
    pub fn j_label_rule(&self) -> &str {
        &self.j_label_rule
    }
    pub fn sources(&self) -> &[Source] {
        &self.sources
    }
    pub fn v_entries(&self) -> &[VReference] {
        &self.v_entries
    }
    pub fn j_entries(&self) -> &[JReference] {
        &self.j_entries
    }
    pub fn revision(&self) -> Result<String> {
        self.validate()?;
        hash(
            b"antibody-reference-corpus-3\0",
            &(REFERENCE_ENGINE, self),
            1 << 20,
        )
    }
    fn validate(&self) -> Result<()> {
        // Native producers receive a serialized-byte bound before cloning any
        // labels/ledger facts during semantic validation.
        struct Count(usize);
        impl std::io::Write for Count {
            fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
                self.0 = self
                    .0
                    .checked_sub(bytes.len())
                    .ok_or_else(|| std::io::Error::other("reference corpus bytes"))?;
                Ok(bytes.len())
            }
            fn flush(&mut self) -> std::io::Result<()> {
                Ok(())
            }
        }
        serde_json::to_writer(Count(1 << 20), self)
            .map_err(|_| Error::Limit("embedded reference bytes"))?;
        if self.version != 3
            || self.sources.len() > 3
            || self.v_entries.len() > 354
            || self.j_entries.len() > 23
        {
            return Err(Error::Invalid("bounded current reference inventory"));
        }
        let (scope, expected, v_total, j_total) = match (&self.species, &self.source_policy) {
            (Species::Human, SourcePolicy::HumanProviderLabeled) => (
                "human-functional-V-prefix-and-J-FR4",
                vec![
                    ("OGRDB:G00060", 10, "IGH", 194, 7),
                    ("OGRDB:G00061", 5, "IGK", 74, 7),
                    ("OGRDB:G00062", 4, "IGL", 86, 9),
                ],
                354,
                23,
            ),
            (
                Species::MouseC57Bl6,
                SourcePolicy::MousePinnedColumns {
                    gapped_sha256,
                    v_label_evidence,
                    atypical_anchors,
                    j_phase_policy,
                },
            ) => {
                if gapped_sha256
                    != "27403baa82141e88a065ddedbc9fc6b4bc4594fc9ecf64d8c4056e98faf4d160"
                    || v_label_evidence
                        != "Ordinals 1..104 of the paired provider-declared IMGT-gapped triplet columns, crosschecked against exact source coding DNA and CDR nucleotide intervals; JSON labels are absent."
                    || j_phase_policy
                        != "Pinned-source observation: j_cdr3_end-1 locates W and has phase frame-1. The AIRR textual frame interpretation (4-frame)%3 differs for these records; both are retained. No generic frame reinterpretation or J numbering labels are supplied."
                    || atypical_anchors.len() != 4
                {
                    return Err(Error::Invalid("mouse column/anchor/frame provenance"));
                }
                let mut observed = Vec::new();
                for entry in &self.v_entries {
                    if entry.aligned_protein.len() != 104
                        || entry.aligned_nucleotides.len() != 312
                        || !entry.aligned_protein.is_ascii()
                        || !entry.aligned_nucleotides.is_ascii()
                    {
                        return Err(Error::Invalid("mouse reference columns"));
                    }
                    for (column, canonical) in [(23_u16, 'C'), (41, 'W'), (104, 'C')] {
                        let index = usize::from(column) - 1;
                        let residue = entry.aligned_protein.as_bytes()[index] as char;
                        if residue != canonical {
                            observed.push(AtypicalAnchor {
                                allele_ref: entry.allele_ref.clone(),
                                imgt_column: column,
                                codon: entry.aligned_nucleotides[index * 3..index * 3 + 3].into(),
                                residue,
                                canonical_residue: canonical,
                            });
                        }
                    }
                }
                observed.sort_by(|a, b| {
                    (&a.allele_ref, a.imgt_column).cmp(&(&b.allele_ref, b.imgt_column))
                });
                if &observed != atypical_anchors {
                    return Err(Error::Invalid("mouse functional atypical anchor ledger"));
                }
                (
                    "mouse-C57BL6-functional-IGH-V-prefix-and-observed-J-FR4",
                    vec![("OGRDB:G00004", 5, "IGH", 95, 4)],
                    95,
                    4,
                )
            }
            _ => return Err(Error::Invalid("species/reference-policy association")),
        };
        if self.version != 3
            || self.license != "CC0-1.0"
            || self.license_url
                != "https://ogrdb.airr-community.org/render_page/licensing_statement.html"
            || self.scope != scope
            || self.j_label_rule != "IMGT-118-sequential-complete-J-codons"
            || self.sources.len() != expected.len()
            || self.v_entries.len() != v_total
            || self.j_entries.len() != j_total
        {
            return Err(Error::Invalid("unsupported embedded reference corpus"));
        }
        let mut identities = BTreeSet::new();
        let mut counts = vec![0; expected.len()];
        for (source, &(id, version, locus, retained, j_retained)) in
            self.sources.iter().zip(&expected)
        {
            if source.set_id != id
                || source.release_version != version
                || source.locus != locus
                || source.v_records != retained
                || source.j_records != j_retained
                || source.species
                    != if self.species == Species::Human {
                        "Homo sapiens"
                    } else {
                        "Mus musculus"
                    }
                || source.taxonomy_id
                    != if self.species == Species::Human {
                        "NCBITAXON:9606"
                    } else {
                        "NCBITAXON:10090"
                    }
                || self.species == Species::MouseC57Bl6
                    && (source.sha256
                        != "a557dc47e4946b015fcb5414c29c275783c464848fd13b8f4a0386c4d6e686e8"
                        || source.set_ref != "OGRDB:C57BL/6 IGH.5")
                || source.author.trim().is_empty()
                || source.sha256.len() != 64
                || !source
                    .sha256
                    .bytes()
                    .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
                || source.v_records + source.j_records + source.excluded.len()
                    != source.original_records
            {
                return Err(Error::Invalid("reference source metadata"));
            }
            for excluded in &source.excluded {
                if !(if self.species == Species::Human {
                    excluded.reason == "not-functional-V-or-J"
                } else {
                    matches!(
                        excluded.reason.as_str(),
                        "provider-not-functional" | "D-not-a-V-or-J-reference-model"
                    )
                }) || !identities.insert(excluded.allele_ref.as_str())
                {
                    return Err(Error::Invalid("reference exclusion ledger"));
                }
            }
        }
        for entry in &self.v_entries {
            if entry.source >= self.sources.len()
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
        let mut j_counts = vec![0; expected.len()];
        for entry in &self.j_entries {
            if entry.source >= self.sources.len()
                || entry.release_version == 0
                || entry.label.is_empty()
                || entry.allele_id.is_empty()
                || !identities.insert(entry.allele_ref.as_str())
                || entry.original_sequence.len() > 512
                || entry.coding_nucleotides.len() > 256
                || entry.coding_nucleotides.is_empty()
                || entry.gene_start == 0
                || entry
                    .original_sequence
                    .get(entry.gene_start - 1..entry.gene_end)
                    != Some(entry.coding_nucleotides.as_str())
                || !entry
                    .original_sequence
                    .bytes()
                    .all(|b| b"ACGT".contains(&b))
                || !(1..=3).contains(&entry.provider_codon_frame)
                || entry.anchor_offset >= entry.coding_nucleotides.len()
                || entry.provider_j_cdr3_end != entry.anchor_offset + 1
                || entry.anchor_offset % 3 != usize::from(entry.provider_codon_frame - 1)
            {
                return Err(Error::Invalid("J source interval/anchor/frame"));
            }
            let expected_length = if entry.source == 0 { 11 } else { 10 };
            let end = entry.anchor_offset + entry.fr4_nucleotides.len();
            if entry.fr4_protein.len() != expected_length
                || entry.fr4_nucleotides.len() != expected_length * 3
                || entry.fr4_positions
                    != (118..118 + expected_length)
                        .map(|n| n.to_string())
                        .collect::<Vec<_>>()
                || entry.coding_nucleotides.get(entry.anchor_offset..end)
                    != Some(entry.fr4_nucleotides.as_str())
                || entry.coding_nucleotides.get(end..) != Some(entry.trailing_nucleotides.as_str())
                || entry.trailing_nucleotides.len() != 1
                || !entry
                    .fr4_protein
                    .bytes()
                    .all(|b| PROTEIN_RESIDUES.contains(&b))
            {
                return Err(Error::Invalid("J FR4 complete-codon extent/labels"));
            }
            let aa = entry.fr4_protein.as_bytes();
            if aa[0] != if entry.source == 0 { b'W' } else { b'F' }
                || aa[1] != b'G'
                || aa[3] != b'G'
            {
                return Err(Error::Invalid("J FR4 anchor motif"));
            }
            j_counts[entry.source] += 1;
        }
        if j_counts != expected.iter().map(|e| e.4).collect::<Vec<_>>() {
            return Err(Error::Invalid("J reference source counts"));
        }
        if counts != expected.iter().map(|e| e.3).collect::<Vec<_>>() {
            return Err(Error::Invalid("reference source counts"));
        }
        Ok(())
    }
}
/// Load a finite species-scoped compact corpus once. Raw OGRDB parsing is
/// maintainer/test-only, never a dropdown or application analysis operation.
pub fn load(species: &Species) -> Result<&'static Corpus> {
    static HUMAN: std::sync::OnceLock<Corpus> = std::sync::OnceLock::new();
    static MOUSE: std::sync::OnceLock<Corpus> = std::sync::OnceLock::new();
    let (slot, data) = match species {
        Species::Human => (&HUMAN, HUMAN_DATA),
        Species::MouseC57Bl6 => (&MOUSE, MOUSE_DATA),
    };
    if let Some(corpus) = slot.get() {
        return Ok(corpus);
    }
    if data.len() > 1 << 20 {
        return Err(Error::Limit("embedded reference bytes"));
    }
    let corpus: Corpus =
        serde_json::from_str(data).map_err(|_| Error::Invalid("embedded reference JSON"))?;
    corpus.validate()?;
    if &corpus.species != species {
        return Err(Error::Binding);
    }
    let _ = slot.set(corpus);
    slot.get().ok_or(Error::Binding)
}
