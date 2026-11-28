//! Observed reference-CDS usage. Counts/preferences are not expression calibration.
//! Derivation is curation work; the runtime catalogue embeds the compact profile
//! and never parses a genome for a dropdown interaction.
mod input;
mod runtime;
pub use runtime::{Reference, Summary, catalogue, resolve};
#[cfg(test)]
mod tests;

use super::CodonWeight;
use crate::sequence::translation::{self, CodeDescriptor, InitiationPolicy, TranslationOptions};
use collagene_model::{Error, Sequence};
use serde::{Deserialize, Serialize};
use std::io::Read;

pub const ENGINE: &str = "reference-codon-counts-1";
pub const POLICY: &str = "complete-chromosome-internal-cds-1";
pub const VERSION: u16 = 1;
pub const MAX_COMPRESSED_BYTES: usize = 20 * 1024 * 1024;
pub const MAX_EXPANDED_BYTES: usize = 32 * 1024 * 1024;
pub const MAX_RECORDS: usize = 20_000;
pub const MAX_CDS_BASES: usize = 300_000;
pub const MAX_LINE_BYTES: usize = 4096;
const URL_PREFIX: &str =
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/";

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Dataset {
    pub id: String,
    pub organism: String,
    pub assembly: String,
    pub assembly_name: String,
    pub taxid: u32,
    pub sequence: String,
    pub sequence_length: usize,
    pub genetic_code: u8,
    pub retrieved_utc_date: String,
    pub source_url_prefix: String,
}
impl Dataset {
    pub fn mg1655() -> Self {
        Self {
            id: "ecoli-k12-mg1655-reference".into(),
            organism: "Escherichia coli str. K-12 substr. MG1655".into(),
            assembly: "GCF_000005845.2".into(),
            assembly_name: "ASM584v2".into(),
            taxid: 511145,
            sequence: "NC_000913.3".into(),
            sequence_length: 4_641_652,
            genetic_code: 11,
            retrieved_utc_date: "2026-11-28".into(),
            source_url_prefix: URL_PREFIX.into(),
        }
    }
    fn validate(&self) -> Result<(), Error> {
        for value in [
            &self.id,
            &self.organism,
            &self.assembly,
            &self.assembly_name,
            &self.sequence,
            &self.retrieved_utc_date,
            &self.source_url_prefix,
        ] {
            if value.is_empty() || value.len() > 512 || value.chars().any(char::is_control) {
                return Err(bad("invalid dataset identity/provenance text"));
            }
        }
        if self.genetic_code != 11
            || self.taxid == 0
            || self.sequence_length == 0
            || self.sequence_length > collagene_model::MAX_SEQUENCE_LEN
            || !self.source_url_prefix.starts_with("https://")
            || !self.source_url_prefix.ends_with('/')
        {
            return Err(bad(
                "this derivation policy requires a bounded chromosome and genetic code 11",
            ));
        }
        CodeDescriptor::for_id(self.genetic_code)?;
        Ok(())
    }
}

/// Raw byte identity, distinct from the annotation's mutable release label.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Source {
    pub url: String,
    pub compressed: bool,
    pub bytes: u64,
    pub blake3: String,
    pub expanded_bytes: u64,
    pub expanded_blake3: String,
}
impl Source {
    fn new(dataset: &Dataset, suffix: &str, raw: &[u8], expanded: &[u8], compressed: bool) -> Self {
        Self {
            url: format!("{}{}", dataset.source_url_prefix, suffix),
            compressed,
            bytes: raw.len() as u64,
            blake3: hash(raw),
            expanded_bytes: expanded.len() as u64,
            expanded_blake3: hash(expanded),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum Exclusion {
    Pseudogene,
    TranslationException,
    Partial,
    OtherGeneticCode,
    MissingProtein,
    NonzeroPhase,
    Multipart,
    AmbiguousBases,
    PartialCodon,
    UnrecognizedStart,
    MissingTerminalStop,
    InternalStop,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Excluded {
    pub reason: Exclusion,
    pub loci: u64,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct HostProfile {
    pub version: u16,
    pub engine: String,
    pub policy: String,
    pub dataset: Dataset,
    /// CDS FASTA, genomic GFF3, assembly report, annotation hashes, in that order.
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Source, 4>")]
    pub sources: Vec<Source>,
    /// Exact retained metadata text, including annotation-change timestamps.
    pub annotation_hashes: String,
    /// A,C,G,T lexical codon order; only internal ordinary sense codons counted.
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, u64, 64>")]
    pub counts: Vec<u64>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, u64, 64>")]
    pub initiation_counts: Vec<u64>,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, u64, 64>")]
    pub terminal_stop_counts: Vec<u64>,
    pub observed_loci: u64,
    pub admitted_loci: u64,
    pub internal_codons: u64,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Excluded, 12>")]
    pub exclusions: Vec<Excluded>,
    /// Sorted admitted chromosome/locus/protein/coordinate identities, not proteins.
    pub admitted_loci_blake3: String,
    pub identity: String,
}

#[derive(Clone, Copy)]
pub struct Expanded<'a> {
    pub cds: &'a [u8],
    pub gff: &'a [u8],
    pub assembly_report: &'a [u8],
    pub annotation_hashes: &'a [u8],
}
#[derive(Clone, Copy)]
pub struct Compressed<'a> {
    pub cds_gzip: &'a [u8],
    pub gff_gzip: &'a [u8],
    pub assembly_report: &'a [u8],
    pub annotation_hashes: &'a [u8],
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Usage {
    pub codon: String,
    pub amino_acid: u8,
    pub count: u64,
    /// Integer millionths of all observed internal sense codons.
    pub total_millionths: u32,
    /// None means no observation in this synonymous family, not a forbidden codon.
    pub synonymous_millionths: Option<u32>,
    pub preference: u32,
}
impl HostProfile {
    /// Fresh curation/test derivation. Not the eventual runtime catalogue accessor.
    pub fn derive_mg1655() -> Result<Self, Error> {
        Self::derive_gzip(
            Dataset::mg1655(),
            Compressed {
                cds_gzip: include_bytes!("data/mg1655/cds.fna.gz"),
                gff_gzip: include_bytes!("data/mg1655/genomic.gff.gz"),
                assembly_report: include_bytes!("data/mg1655/assembly_report.txt"),
                annotation_hashes: include_bytes!("data/mg1655/annotation_hashes.txt"),
            },
        )
    }
    pub fn derive_gzip(dataset: Dataset, raw: Compressed<'_>) -> Result<Self, Error> {
        dataset.validate()?;
        if raw
            .cds_gzip
            .len()
            .checked_add(raw.gff_gzip.len())
            .is_none_or(|n| n > MAX_COMPRESSED_BYTES)
        {
            return Err(Error::TooLarge);
        }
        check_metadata(raw.assembly_report, raw.annotation_hashes)?;
        let cds = expand(raw.cds_gzip, MAX_EXPANDED_BYTES)?;
        let gff = expand(raw.gff_gzip, MAX_EXPANDED_BYTES.saturating_sub(cds.len()))?;
        let sources = source_records(
            &dataset,
            raw.cds_gzip,
            raw.gff_gzip,
            &cds,
            &gff,
            raw.assembly_report,
            raw.annotation_hashes,
            true,
        );
        Self::derive(
            dataset,
            Expanded {
                cds: &cds,
                gff: &gff,
                assembly_report: raw.assembly_report,
                annotation_hashes: raw.annotation_hashes,
            },
            sources,
        )
    }
    pub fn derive_expanded(dataset: Dataset, raw: Expanded<'_>) -> Result<Self, Error> {
        dataset.validate()?;
        check_metadata(raw.assembly_report, raw.annotation_hashes)?;
        if raw
            .cds
            .len()
            .checked_add(raw.gff.len())
            .is_none_or(|n| n > MAX_EXPANDED_BYTES)
        {
            return Err(Error::TooLarge);
        }
        let sources = source_records(
            &dataset,
            raw.cds,
            raw.gff,
            raw.cds,
            raw.gff,
            raw.assembly_report,
            raw.annotation_hashes,
            false,
        );
        Self::derive(dataset, raw, sources)
    }
    fn derive(dataset: Dataset, raw: Expanded<'_>, sources: Vec<Source>) -> Result<Self, Error> {
        check_metadata(raw.assembly_report, raw.annotation_hashes)?;
        if raw
            .cds
            .len()
            .checked_add(raw.gff.len())
            .is_none_or(|n| n > MAX_EXPANDED_BYTES)
        {
            return Err(Error::TooLarge);
        }
        input::validate_metadata(&dataset, raw.assembly_report, raw.annotation_hashes)?;
        let gff = input::Annotation::read(raw.gff, &dataset)?;
        let records = input::CodingSequences::read(raw.cds, &dataset, &gff)?;
        let mut profile = Self {
            version: VERSION,
            engine: ENGINE.into(),
            policy: POLICY.into(),
            dataset,
            sources,
            annotation_hashes: text(raw.annotation_hashes)?.into(),
            counts: vec![0; 64],
            initiation_counts: vec![0; 64],
            terminal_stop_counts: vec![0; 64],
            observed_loci: records.rows.len() as u64,
            admitted_loci: 0,
            internal_codons: 0,
            exclusions: vec![],
            admitted_loci_blake3: String::new(),
            identity: String::new(),
        };
        let mut excluded = std::collections::BTreeMap::<Exclusion, u64>::new();
        let mut admitted = blake3::Hasher::new();
        let code = CodeDescriptor::for_id(profile.dataset.genetic_code)?;
        for record in &records.rows {
            if let Some(reason) = record.exclusion(&code)? {
                *excluded.entry(reason).or_default() += 1;
                continue;
            }
            profile.admitted_loci += 1;
            profile.initiation_counts[index(&record.bases.as_bytes()[..3])?] += 1;
            profile.terminal_stop_counts
                [index(&record.bases.as_bytes()[record.bases.len() - 3..])?] += 1;
            for codon in record.bases.as_bytes()[3..record.bases.len() - 3].chunks_exact(3) {
                profile.counts[index(codon)?] += 1;
                profile.internal_codons += 1;
            }
            // Input owner orders rows by the full locus identity, independent of file row order.
            admitted.update(&(record.identity.len() as u64).to_le_bytes());
            admitted.update(record.identity.as_bytes());
        }
        profile.exclusions = excluded
            .into_iter()
            .map(|(reason, loci)| Excluded { reason, loci })
            .collect();
        profile.admitted_loci_blake3 = admitted.finalize().to_hex().to_string();
        profile.identity = profile.expected_identity()?;
        profile.validate()?;
        Ok(profile)
    }
    fn expected_identity(&self) -> Result<String, Error> {
        // Explicit fields exclude only the identity itself; no second representation.
        let bytes = serde_json::to_vec(&(
            self.version,
            &self.engine,
            &self.policy,
            &self.dataset,
            &self.sources,
            &self.annotation_hashes,
            &self.counts,
            &self.initiation_counts,
            &self.terminal_stop_counts,
            self.observed_loci,
            self.admitted_loci,
            self.internal_codons,
            &self.exclusions,
            &self.admitted_loci_blake3,
        ))
        .map_err(|_| Error::TooLarge)?;
        if bytes.len() > 32 * 1024 {
            return Err(Error::TooLarge);
        }
        Ok(hash(&bytes))
    }
    pub fn validate(&self) -> Result<(), Error> {
        if self.version != VERSION || self.engine != ENGINE || self.policy != POLICY {
            return Err(Error::Unsupported(
                "host codon profile/policy version".into(),
            ));
        }
        self.dataset.validate()?;
        if self.sources.len() != 4
            || self.annotation_hashes.len() > 8192
            || self.counts.len() != 64
            || self.initiation_counts.len() != 64
            || self.terminal_stop_counts.len() != 64
            || self.observed_loci == 0
            || self.observed_loci > MAX_RECORDS as u64
            || self.admitted_loci == 0
            || self.admitted_loci > self.observed_loci
            || self.internal_codons == 0
            || self.internal_codons > MAX_EXPANDED_BYTES as u64 / 3
            || !digest(&self.admitted_loci_blake3)
            || !digest(&self.identity)
        {
            return Err(bad("invalid profile dimensions/counts"));
        }
        let base = format!("{}_{}", self.dataset.assembly, self.dataset.assembly_name);
        let compressed = self.sources[0].compressed;
        let suffix = if compressed { ".gz" } else { "" };
        let expected = [
            format!("{base}_cds_from_genomic.fna{suffix}"),
            format!("{base}_genomic.gff{suffix}"),
            format!("{base}_assembly_report.txt"),
            "annotation_hashes.txt".into(),
        ];
        for (index, source) in self.sources.iter().enumerate() {
            if source.url != format!("{}{}", self.dataset.source_url_prefix, expected[index])
                || !digest(&source.blake3)
                || !digest(&source.expanded_blake3)
                || source.bytes == 0
                || source.expanded_bytes == 0
                || source.bytes > MAX_EXPANDED_BYTES as u64
                || source.expanded_bytes > MAX_EXPANDED_BYTES as u64
                || source.compressed != (index < 2 && compressed)
                || (!source.compressed
                    && (source.bytes != source.expanded_bytes
                        || source.blake3 != source.expanded_blake3))
            {
                return Err(bad("invalid source provenance"));
            }
        }
        if self.sources[0].expanded_bytes + self.sources[1].expanded_bytes
            > MAX_EXPANDED_BYTES as u64
            || (compressed
                && self.sources[0].bytes + self.sources[1].bytes > MAX_COMPRESSED_BYTES as u64)
            || self.sources[2].bytes > 16 * 1024
            || self.sources[3].bytes > 8192
            || self.sources[3].blake3 != hash(self.annotation_hashes.as_bytes())
        {
            return Err(bad("source sizes/retained annotation hash disagree"));
        }
        let sum = |values: &[u64]| {
            values.iter().try_fold(0u64, |sum, value| {
                sum.checked_add(*value).ok_or(Error::TooLarge)
            })
        };
        if sum(&self.counts)? != self.internal_codons
            || sum(&self.initiation_counts)? != self.admitted_loci
            || sum(&self.terminal_stop_counts)? != self.admitted_loci
        {
            return Err(bad("profile count sums disagree"));
        }
        let mut previous = None;
        let mut rejected = 0u64;
        for exclusion in &self.exclusions {
            if exclusion.loci == 0 || previous.is_some_and(|p| p >= exclusion.reason) {
                return Err(bad("invalid exclusion order/count"));
            }
            previous = Some(exclusion.reason);
            rejected = rejected
                .checked_add(exclusion.loci)
                .ok_or(Error::TooLarge)?;
        }
        if self.admitted_loci.checked_add(rejected) != Some(self.observed_loci) {
            return Err(bad("unaccounted CDS records"));
        }
        let amino = amino_acids(self.dataset.genetic_code)?;
        let starts = CodeDescriptor::for_id(self.dataset.genetic_code)?.start_codons;
        for i in 0..64 {
            if (amino[i] == b'*' && self.counts[i] != 0)
                || (amino[i] != b'*' && self.terminal_stop_counts[i] != 0)
                || (self.initiation_counts[i] != 0 && !starts.iter().any(|s| s == &codon(i)))
            {
                return Err(bad("counts disagree with canonical genetic code"));
            }
        }
        if self.expected_identity()? != self.identity {
            return Err(bad("profile identity differs"));
        }
        Ok(())
    }
    pub fn usage(&self) -> Result<Vec<Usage>, Error> {
        self.validate()?;
        let amino = amino_acids(self.dataset.genetic_code)?;
        let mut rows = Vec::with_capacity(64);
        for i in 0..64 {
            let sum: u64 = (0..64)
                .filter(|j| amino[*j] == amino[i])
                .map(|j| self.counts[j])
                .sum();
            let maximum = (0..64)
                .filter(|j| amino[*j] == amino[i])
                .map(|j| self.counts[j])
                .max()
                .unwrap_or(0);
            let count = self.counts[i];
            rows.push(Usage {
                codon: codon(i),
                amino_acid: amino[i],
                count,
                total_millionths: (count * 1_000_000 / self.internal_codons) as u32,
                synonymous_millionths: if sum == 0 {
                    None
                } else {
                    Some((count * 1_000_000 / sum) as u32)
                },
                preference: if maximum == 0 || amino[i] == b'*' {
                    1
                } else {
                    1 + (count * 999_999 / maximum) as u32
                },
            });
        }
        Ok(rows)
    }
    /// All 64 explicit nonzero weights. Unobserved codons remain allowed.
    /// Arithmetic score preferences, not CAI or expression/synthesis probability.
    pub fn preferences(&self) -> Result<Vec<CodonWeight>, Error> {
        Ok(self
            .usage()?
            .into_iter()
            .map(|row| CodonWeight {
                codon: row.codon,
                weight: row.preference,
            })
            .collect())
    }
}
fn source_records(
    dataset: &Dataset,
    raw_cds: &[u8],
    raw_gff: &[u8],
    cds: &[u8],
    gff: &[u8],
    report: &[u8],
    annotations: &[u8],
    compressed: bool,
) -> Vec<Source> {
    let base = format!("{}_{}", dataset.assembly, dataset.assembly_name);
    vec![
        Source::new(
            dataset,
            &format!(
                "{base}_cds_from_genomic.fna{}",
                if compressed { ".gz" } else { "" }
            ),
            raw_cds,
            cds,
            compressed,
        ),
        Source::new(
            dataset,
            &format!("{base}_genomic.gff{}", if compressed { ".gz" } else { "" }),
            raw_gff,
            gff,
            compressed,
        ),
        Source::new(
            dataset,
            &format!("{base}_assembly_report.txt"),
            report,
            report,
            false,
        ),
        Source::new(
            dataset,
            "annotation_hashes.txt",
            annotations,
            annotations,
            false,
        ),
    ]
}
fn expand(bytes: &[u8], maximum: usize) -> Result<Vec<u8>, Error> {
    // MultiGzDecoder checks all gzip members and trailing junk; the read cap bounds
    // expansion before any input-controlled reserve or allocation.
    let mut out = Vec::new();
    flate2::read::MultiGzDecoder::new(bytes)
        .take(maximum as u64 + 1)
        .read_to_end(&mut out)
        .map_err(bad)?;
    if out.len() > maximum {
        return Err(Error::TooLarge);
    }
    Ok(out)
}
fn check_metadata(report: &[u8], annotations: &[u8]) -> Result<(), Error> {
    if report.len() > 16 * 1024 || annotations.len() > 8192 {
        Err(Error::TooLarge)
    } else {
        Ok(())
    }
}
fn amino_acids(code: u8) -> Result<Vec<u8>, Error> {
    let all = (0..64).map(codon).collect::<String>();
    Ok(translation::translate(
        &Sequence::dna(&all)?,
        &TranslationOptions {
            genetic_code: code,
            ..Default::default()
        },
    )?
    .protein
    .into_bytes())
}
fn codon(index: usize) -> String {
    [
        b"ACGT"[index / 16],
        b"ACGT"[(index / 4) % 4],
        b"ACGT"[index % 4],
    ]
    .into_iter()
    .map(char::from)
    .collect()
}
fn index(codon: &[u8]) -> Result<usize, Error> {
    codon.iter().try_fold(0usize, |index, base| {
        Ok(index * 4
            + b"ACGT"
                .iter()
                .position(|b| b == base)
                .ok_or_else(|| bad("nonconcrete counted codon"))?)
    })
}
fn hash(bytes: &[u8]) -> String {
    blake3::hash(bytes).to_hex().to_string()
}
fn digest(text: &str) -> bool {
    text.len() == 64
        && text
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn text(bytes: &[u8]) -> Result<&str, Error> {
    std::str::from_utf8(bytes).map_err(bad)
}
fn bad(error: impl std::fmt::Display) -> Error {
    Error::Format(format!("host codon profile: {error}"))
}
