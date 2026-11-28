//! Pinned chromosome CDS↔GFF locus binding; no genomic-sequence or expression claim.
use super::*;
use collagene_model::Location;
use std::collections::BTreeMap;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Key {
    sequence: String,
    locus: String,
    protein: String,
}
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct Part {
    start: usize,
    end: usize,
    reverse: bool,
}
#[derive(Debug, Clone)]
struct Gene {
    sequence: String,
    locus: String,
    start: usize,
    end: usize,
    reverse: bool,
    pseudo: bool,
    partial: bool,
}
struct CodingFeature {
    id: String,
    parent: String,
    parts: Vec<Part>,
    pseudo: bool,
    partial: bool,
    exception: bool,
    code: Option<u8>,
    phase_zero: bool,
}
pub(super) struct Annotation {
    groups: BTreeMap<Key, CodingFeature>,
}
impl Annotation {
    pub fn read(bytes: &[u8], dataset: &Dataset) -> Result<Self, Error> {
        let source = text(bytes)?;
        let mut groups = BTreeMap::<Key, CodingFeature>::new();
        let mut genes = BTreeMap::<String, Gene>::new();
        let mut ids = BTreeMap::<String, Key>::new();
        let mut rows = 0usize;
        let mut assembly = None;
        let mut sequence_region = false;
        for line in source.lines() {
            if line.len() > MAX_LINE_BYTES {
                return Err(Error::TooLarge);
            }
            if let Some(value) = line.strip_prefix("#!genome-build-accession NCBI_Assembly:") {
                if assembly.replace(value).is_some() {
                    return Err(bad("duplicate GFF assembly directive"));
                }
            }
            if let Some(value) = line.strip_prefix("##sequence-region ") {
                let fields = value.split_whitespace().collect::<Vec<_>>();
                if fields.len() != 3
                    || fields[0] != dataset.sequence
                    || fields[1] != "1"
                    || integer(fields[2])? != dataset.sequence_length
                {
                    return Err(bad("GFF chromosome region differs from declared dataset"));
                }
                if sequence_region {
                    return Err(bad("duplicate GFF sequence region"));
                }
                sequence_region = true;
            }
            if line.starts_with("##FASTA") {
                return Err(bad("embedded GFF FASTA is outside this CDS-source policy"));
            }
            if line.starts_with('#') || line.is_empty() {
                continue;
            }
            rows += 1;
            if rows > 100_000 {
                return Err(Error::TooLarge);
            }
            let fields = line.split('\t').collect::<Vec<_>>();
            if fields.len() != 9 {
                return Err(bad("GFF requires nine fields"));
            }
            if !matches!(fields[2], "gene" | "CDS" | "pseudogene") {
                continue;
            }
            if fields[0] != dataset.sequence {
                return Err(bad("GFF feature belongs to another chromosome"));
            }
            let start = integer(fields[3])?
                .checked_sub(1)
                .ok_or_else(|| bad("GFF coordinate zero"))?;
            let end = integer(fields[4])?;
            if start >= end || end > dataset.sequence_length {
                return Err(bad("GFF feature coordinates"));
            }
            let reverse = match fields[6] {
                "+" => false,
                "-" => true,
                _ => return Err(bad("CDS/gene strand must be explicit")),
            };
            let attrs = attributes(fields[8])?;
            let id = required(&attrs, "ID")?;
            let locus = required(&attrs, "locus_tag")?;
            identifier(id)?;
            identifier(locus)?;
            let pseudo = boolean(&attrs, "pseudo")?
                || attrs.contains_key("pseudogene")
                || attrs.get("gene_biotype") == Some(&"pseudogene")
                || fields[2] == "pseudogene";
            let partial = boolean(&attrs, "partial")?
                || attrs.contains_key("start_range")
                || attrs.contains_key("end_range");
            if fields[2] != "CDS" {
                if genes.len() >= MAX_RECORDS
                    || genes
                        .insert(
                            id.into(),
                            Gene {
                                sequence: fields[0].into(),
                                locus: locus.into(),
                                start,
                                end,
                                reverse,
                                pseudo,
                                partial,
                            },
                        )
                        .is_some()
                {
                    return Err(bad("duplicate/excessive GFF gene identities"));
                }
                continue;
            }
            let parent = required(&attrs, "Parent")?;
            identifier(parent)?;
            let protein = attrs.get("protein_id").copied().unwrap_or("");
            if !protein.is_empty() {
                identifier(protein)?;
            }
            let code = attrs
                .get("transl_table")
                .map(|v| v.parse::<u8>().map_err(bad))
                .transpose()?;
            let phase_zero = match fields[7] {
                "0" => true,
                "1" | "2" | "." => false,
                _ => return Err(bad("invalid GFF CDS phase")),
            };
            let exception = attrs.contains_key("exception") || attrs.contains_key("transl_except");
            let key = Key {
                sequence: fields[0].into(),
                locus: locus.into(),
                protein: protein.into(),
            };
            if let Some(previous) = ids.insert(id.into(), key.clone()) {
                if previous != key {
                    return Err(bad("GFF CDS ID reused by another locus/protein"));
                }
            }
            if !groups.contains_key(&key) && groups.len() >= MAX_RECORDS {
                return Err(Error::TooLarge);
            }
            let group = groups.entry(key).or_insert_with(|| CodingFeature {
                id: id.into(),
                parent: parent.into(),
                parts: vec![],
                pseudo,
                partial,
                exception,
                code,
                phase_zero,
            });
            if group.id != id
                || group.parent != parent
                || group.code != code
                || group.pseudo != pseudo
                || group.partial != partial
                || group.exception != exception
            {
                return Err(bad("inconsistent grouped GFF CDS identity/flags"));
            }
            if group.parts.len() >= 64 {
                return Err(Error::TooLarge);
            }
            group.phase_zero &= phase_zero;
            group.parts.push(Part {
                start,
                end,
                reverse,
            });
        }
        if assembly != Some(dataset.assembly.as_str()) || !sequence_region || groups.is_empty() {
            return Err(bad(
                "GFF assembly/sequence declaration is missing or differs",
            ));
        }
        for (key, group) in &mut groups {
            let gene = genes
                .get(&group.parent)
                .ok_or_else(|| bad("CDS parent gene unavailable"))?;
            if gene.sequence != key.sequence
                || gene.locus != key.locus
                || group
                    .parts
                    .iter()
                    .any(|p| p.start < gene.start || p.end > gene.end || p.reverse != gene.reverse)
            {
                return Err(bad("CDS locus/coordinates differ from parent gene"));
            }
            group.pseudo |= gene.pseudo;
            group.partial |= gene.partial;
            group.parts.sort();
            if group.parts.windows(2).any(|p| p[0] == p[1]) {
                return Err(bad("duplicate GFF CDS coordinate row"));
            }
        }
        Ok(Self { groups })
    }
}

pub(super) struct CodingSequences {
    pub rows: Vec<CodingSequence>,
}
pub(super) struct CodingSequence {
    pub identity: String,
    pub bases: String,
    pseudo: bool,
    exception: bool,
    partial: bool,
    code: Option<u8>,
    phase_zero: bool,
    protein_present: bool,
    multipart: bool,
}
impl CodingSequences {
    pub fn read(bytes: &[u8], dataset: &Dataset, gff: &Annotation) -> Result<Self, Error> {
        let source = text(bytes)?;
        let mut records = BTreeMap::<Key, CodingSequence>::new();
        let mut header = None;
        let mut bases = String::new();
        for line in source.lines() {
            if line.len() > MAX_LINE_BYTES {
                return Err(Error::TooLarge);
            }
            if line.starts_with('>') {
                if let Some(previous) = header.replace(line) {
                    Self::record(
                        previous,
                        std::mem::take(&mut bases),
                        dataset,
                        gff,
                        &mut records,
                    )?;
                }
            } else {
                if header.is_none()
                    || line.is_empty()
                    || !line.bytes().all(|b| b"ACGTRYSWKMBDHVN".contains(&b))
                {
                    return Err(bad("invalid CDS FASTA sequence line"));
                }
                if bases.len() + line.len() > MAX_CDS_BASES {
                    return Err(Error::TooLarge);
                }
                bases.push_str(line);
            }
        }
        if let Some(previous) = header {
            Self::record(previous, bases, dataset, gff, &mut records)?;
        }
        if records.len() != gff.groups.len() || records.is_empty() {
            return Err(bad("CDS FASTA and GFF loci are not one-to-one"));
        }
        Ok(Self {
            rows: records.into_values().collect(),
        })
    }
    fn record(
        header: &str,
        bases: String,
        dataset: &Dataset,
        gff: &Annotation,
        records: &mut BTreeMap<Key, CodingSequence>,
    ) -> Result<(), Error> {
        if records.len() >= MAX_RECORDS {
            return Err(Error::TooLarge);
        }
        let (id, _) = header
            .split_once(' ')
            .ok_or_else(|| bad("CDS header attributes missing"))?;
        let (sequence, local_identity) = id
            .strip_prefix(">lcl|")
            .and_then(|s| s.split_once("_cds_"))
            .ok_or_else(|| bad("CDS local sequence identity"))?;
        if sequence != dataset.sequence {
            return Err(bad("CDS header chromosome differs"));
        }
        let attrs = fasta_attributes(header)?;
        if attrs.get("gbkey") != Some(&"CDS") {
            return Err(bad("FASTA record is not declared CDS"));
        }
        let locus = required(&attrs, "locus_tag")?;
        let protein = attrs.get("protein_id").copied().unwrap_or("");
        identifier(locus)?;
        if !protein.is_empty() {
            identifier(protein)?;
        }
        let declared_identity = if protein.is_empty() { locus } else { protein };
        let ordinal = local_identity
            .strip_prefix(&format!("{declared_identity}_"))
            .ok_or_else(|| bad("FASTA local identifier differs from protein/locus attribute"))?;
        if integer(ordinal)? == 0 {
            return Err(bad("FASTA local ordinal must be positive"));
        }
        let key = Key {
            sequence: sequence.into(),
            locus: locus.into(),
            protein: protein.into(),
        };
        let group = gff
            .groups
            .get(&key)
            .ok_or_else(|| bad("CDS locus/protein does not match GFF"))?;
        let location = collagene_model::annotations::insdc::parse_location(
            required(&attrs, "location")?,
            dataset.sequence_length,
        )?;
        let mut parts = vec![];
        let mut partial = false;
        flatten(&location, false, 0, &mut parts, &mut partial)?;
        if parts != group.parts {
            return Err(bad("CDS header and GFF coordinate paths differ"));
        }
        if parts.iter().map(|p| p.end - p.start).sum::<usize>() != bases.len() {
            return Err(bad(
                "CDS sequence length differs from bound coordinate path",
            ));
        }
        if let Some(table) = attrs.get("transl_table") {
            if Some(table.parse::<u8>().map_err(bad)?) != group.code {
                return Err(bad("FASTA/GFF genetic code differs"));
            }
        }
        let pseudo = boolean(&attrs, "pseudo")?;
        let exception = attrs.contains_key("exception") || attrs.contains_key("transl_except");
        // Omission of inherited flags is safe only because the GFF gate also applies.
        if (pseudo && !group.pseudo) || (exception && !group.exception) {
            return Err(bad("FASTA/GFF exclusion flags differ"));
        }
        let identity = serde_json::to_string(&(
            sequence,
            locus,
            protein,
            parts
                .iter()
                .map(|p| (p.start, p.end, p.reverse))
                .collect::<Vec<_>>(),
        ))
        .map_err(bad)?;
        let declared_phase_zero = match attrs.get("codon_start").copied() {
            None | Some("1") => true,
            Some("2" | "3") => false,
            _ => return Err(bad("invalid FASTA codon_start")),
        };
        let record = CodingSequence {
            identity,
            bases,
            pseudo: group.pseudo || pseudo,
            exception: group.exception || exception,
            partial: group.partial || partial || boolean(&attrs, "partial")?,
            code: group.code,
            phase_zero: group.phase_zero && declared_phase_zero,
            protein_present: !protein.is_empty(),
            multipart: parts.len() != 1,
        };
        if records.insert(key, record).is_some() {
            return Err(bad("duplicate CDS locus/protein in FASTA"));
        }
        Ok(())
    }
}
impl CodingSequence {
    pub fn exclusion(&self, code: &CodeDescriptor) -> Result<Option<Exclusion>, Error> {
        let reason = if self.pseudo {
            Some(Exclusion::Pseudogene)
        } else if self.exception {
            Some(Exclusion::TranslationException)
        } else if self.partial {
            Some(Exclusion::Partial)
        } else if self.code != Some(code.id) {
            Some(Exclusion::OtherGeneticCode)
        } else if !self.protein_present {
            Some(Exclusion::MissingProtein)
        } else if !self.phase_zero {
            Some(Exclusion::NonzeroPhase)
        } else if self.multipart {
            Some(Exclusion::Multipart)
        } else if !self.bases.bytes().all(|b| b"ACGT".contains(&b)) {
            Some(Exclusion::AmbiguousBases)
        } else if self.bases.len() < 6 || self.bases.len() % 3 != 0 {
            Some(Exclusion::PartialCodon)
        } else if !code
            .start_codons
            .iter()
            .any(|s| s.as_bytes() == &self.bases.as_bytes()[..3])
        {
            Some(Exclusion::UnrecognizedStart)
        } else {
            None
        };
        if reason.is_some() {
            return Ok(reason);
        }
        let translated = translation::translate(
            &Sequence::dna(&self.bases)?,
            &TranslationOptions {
                genetic_code: code.id,
                frame: 0,
                reverse: false,
                initiation: InitiationPolicy::RequireRecognizedStart,
            },
        )?;
        if !translated.protein.ends_with('*') {
            return Ok(Some(Exclusion::MissingTerminalStop));
        }
        if translated.protein[..translated.protein.len() - 1].contains('*') {
            return Ok(Some(Exclusion::InternalStop));
        }
        if translated.unresolved_codons != 0 || translated.possible_stop_codons != 0 {
            return Err(bad("unexpected unresolved code-11 translation"));
        }
        Ok(None)
    }
}

pub(super) fn validate_metadata(
    dataset: &Dataset,
    report: &[u8],
    annotations: &[u8],
) -> Result<(), Error> {
    let report = text(report)?;
    let field = |prefix: &str| -> Result<&str, Error> {
        let mut found = report
            .lines()
            .filter_map(|line| line.strip_prefix(prefix))
            .map(str::trim);
        let value = found
            .next()
            .ok_or_else(|| bad("assembly report field missing"))?;
        if found.next().is_some() {
            return Err(bad("duplicate assembly report field"));
        }
        Ok(value)
    };
    if field("# Assembly name:")? != dataset.assembly_name
        || field("# Taxid:")?.parse::<u32>().map_err(bad)? != dataset.taxid
        || field("# RefSeq assembly accession:")? != dataset.assembly
        || !field("# Organism name:")?.starts_with(&dataset.organism)
    {
        return Err(bad("assembly report identity differs from dataset"));
    }
    let mut sequences = 0;
    for line in report
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
    {
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() != 10
            || fields[6] != dataset.sequence
            || fields[3] != "Chromosome"
            || integer(fields[8])? != dataset.sequence_length
        {
            return Err(bad("assembly report chromosome identity differs"));
        }
        sequences += 1;
    }
    if sequences != 1 {
        return Err(bad("policy requires one named chromosome"));
    }
    let annotation = text(annotations)?;
    let rows = annotation
        .lines()
        .filter(|line| !line.starts_with('#') && !line.is_empty())
        .collect::<Vec<_>>();
    if rows.len() != 1 || rows[0].split('\t').next() != Some(dataset.assembly.as_str()) {
        return Err(bad("annotation hash record assembly differs"));
    }
    Ok(())
}
fn attributes(text: &str) -> Result<BTreeMap<&str, &str>, Error> {
    let mut result = BTreeMap::new();
    for attribute in text.split(';') {
        let (key, value) = attribute
            .split_once('=')
            .ok_or_else(|| bad("GFF attribute requires key=value"))?;
        if result.len() >= 64
            || key.is_empty()
            || value.is_empty()
            || result.insert(key, value).is_some()
        {
            return Err(bad("duplicate/empty/excessive GFF attributes"));
        }
    }
    Ok(result)
}
fn fasta_attributes(header: &str) -> Result<BTreeMap<&str, &str>, Error> {
    // Product text may itself contain bracketed chemistry such as [4Fe-4S].
    // Only an actual ` [identifier=` delimiter starts the next attribute.
    let mut positions = vec![];
    for (offset, _) in header.match_indices(" [") {
        let rest = &header[offset + 2..];
        let length = rest
            .bytes()
            .take_while(|b| b.is_ascii_alphanumeric() || *b == b'_')
            .count();
        if length != 0 && rest.as_bytes().get(length) == Some(&b'=') {
            if positions.len() >= 64 {
                return Err(Error::TooLarge);
            }
            positions.push((offset, &rest[..length], offset + 2 + length + 1));
        }
    }
    let mut result = BTreeMap::new();
    for (i, (_, key, start)) in positions.iter().enumerate() {
        let end = positions.get(i + 1).map_or(header.len(), |next| next.0);
        let value = header[*start..end]
            .strip_suffix(']')
            .ok_or_else(|| bad("unterminated CDS header attribute"))?;
        if value.is_empty() || result.insert(*key, value).is_some() {
            return Err(bad("duplicate/empty CDS header attribute"));
        }
    }
    Ok(result)
}
fn flatten(
    location: &Location,
    reverse: bool,
    depth: usize,
    parts: &mut Vec<Part>,
    partial: &mut bool,
) -> Result<(), Error> {
    if depth >= 16 || parts.len() >= 64 {
        return Err(Error::TooLarge);
    }
    match location {
        Location::Range {
            start,
            end,
            before,
            after,
        } => {
            *partial |= *before || *after;
            parts.push(Part {
                start: *start,
                end: *end,
                reverse,
            });
        }
        Location::Complement(inner) => flatten(inner, !reverse, depth + 1, parts, partial)?,
        Location::Join(children) => {
            for child in children {
                flatten(child, reverse, depth + 1, parts, partial)?;
            }
        }
        _ => {
            return Err(bad(
                "CDS location is not an exact range/join/complement path",
            ));
        }
    }
    Ok(())
}
fn integer(text: &str) -> Result<usize, Error> {
    if text.is_empty() || !text.bytes().all(|b| b.is_ascii_digit()) {
        return Err(bad("invalid decimal coordinate"));
    }
    text.parse().map_err(bad)
}
fn identifier(text: &str) -> Result<(), Error> {
    if text.is_empty()
        || text.len() > 256
        || !text
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b"_.:|-".contains(&b))
    {
        Err(bad("unsupported locus/protein/parent identity"))
    } else {
        Ok(())
    }
}
fn required<'a>(attrs: &BTreeMap<&str, &'a str>, key: &str) -> Result<&'a str, Error> {
    attrs
        .get(key)
        .copied()
        .ok_or_else(|| bad(format!("required attribute {key} missing")))
}
fn boolean(attrs: &BTreeMap<&str, &str>, key: &str) -> Result<bool, Error> {
    match attrs.get(key).copied() {
        None | Some("false") => Ok(false),
        Some("true") => Ok(true),
        _ => Err(bad("invalid boolean annotation attribute")),
    }
}
