//! GenBank syntax uses the standalone gb-io parser; molecular semantics remain ours.
use crate::model::*;
use base64::{Engine as _, engine::general_purpose::STANDARD};
use gb_io::seq::{self as gb, After, Before};
use serde::{Deserialize, Serialize};
use std::borrow::Cow;

pub const MAX_IMPORT_BYTES: usize = 32 * 1024 * 1024;
const DIGEST_PRODUCER_MARKER: &str =
    "Collagene physical product v1: unverified sequence-only evidence";
const DIGEST_CLAIM_BYTES: usize = 32 * 1024;

fn digest_claim_comments(out: &mut Vec<String>, prefix: &str, value: &str) {
    let mut end = value.len().min(DIGEST_CLAIM_BYTES);
    while !value.is_char_boundary(end) {
        end -= 1;
    }
    let encoded = STANDARD.encode(&value.as_bytes()[..end]);
    let chunks = encoded.len().div_ceil(32);
    for (i, chunk) in encoded.as_bytes().chunks(32).enumerate() {
        out.push(format!(
            "{prefix} {:04}/{chunks:04} {}",
            i + 1,
            std::str::from_utf8(chunk).expect("base64")
        ));
    }
    if end != value.len() {
        out.push("Collagene digest producer data truncated; retain native evidence.".into());
    }
}

/// Exact generated COMMENT/CC marker only; feature notes or prose cannot restore
/// chemical state. Claims stay source data, and the existing unmapped-state guard
/// prevents a sequence-only reimport being silently treated as an intact duplex.
pub(crate) fn retain_digest_producer_comments(m: &mut Molecule, raw: &str) -> bool {
    fn comment(line: &str) -> Option<&str> {
        line.strip_prefix("COMMENT     ")
            .or_else(|| line.strip_prefix("CC   "))
    }
    if !raw
        .lines()
        .filter_map(comment)
        .any(|s| s.trim_end() == DIGEST_PRODUCER_MARKER)
    {
        return false;
    }
    m.metadata
        .insert("source:legacy:chemistry_unmapped".into(), "true".into());
    for (prefix, key) in [
        (
            "CG-DIGEST-RECIPE-v1 ",
            "source:restriction_inspection:recipe",
        ),
        (
            "CG-DIGEST-ENDS-v1 ",
            "source:restriction_inspection:duplex_summary",
        ),
    ] {
        let mut encoded = String::new();
        let mut count = 0usize;
        let mut total = None;
        let mut valid = true;
        for value in raw
            .lines()
            .filter_map(comment)
            .filter_map(|s| s.strip_prefix(prefix))
        {
            count += 1;
            if count > 1366 {
                valid = false;
                break;
            }
            let Some((part, data)) = value.split_once(' ') else {
                valid = false;
                break;
            };
            let Some((index, parts)) = part.split_once('/') else {
                valid = false;
                break;
            };
            let (Ok(index), Ok(parts)) = (index.parse::<usize>(), parts.parse::<usize>()) else {
                valid = false;
                break;
            };
            if index != count
                || parts == 0
                || parts > 1366
                || total.is_some_and(|n| n != parts)
                || data.len() > 32
                || !data
                    .bytes()
                    .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'+' | b'/' | b'='))
            {
                valid = false;
                break;
            }
            total = Some(parts);
            encoded.push_str(data);
        }
        if count == 0 {
            continue;
        }
        let decoded = valid
            .then(|| STANDARD.decode(encoded))
            .transpose()
            .ok()
            .flatten();
        if total == Some(count) {
            if let Some(bytes) = decoded.filter(|b| b.len() <= DIGEST_CLAIM_BYTES) {
                if let Ok(text) = String::from_utf8(bytes) {
                    m.metadata.insert(key.into(), text);
                    continue;
                }
            }
        }
        m.metadata.insert(
            "source:restriction_inspection:invalid_producer_claims".into(),
            "true".into(),
        );
    }
    true
}

// These Collagene primer fields are machine values, not free prose. INSDC physical
// line wrapping must not turn an oligo/hash into different evidence. Only normalize
// a value that validates under that field's declared type; original records remain attached.
fn import_qualifier(key: &str, value: &Option<String>) -> (String, Option<String>) {
    let normalized = value.as_ref().map(|raw| {
        let compact: String = raw.chars().filter(|c| !c.is_ascii_whitespace()).collect();
        match key {
            "source_revision"
            | "library_entry_revision"
            | "library_revision"
            | "library_source_revision"
                if compact.len() == 64 && compact.bytes().all(|c| c.is_ascii_hexdigit()) =>
            {
                compact
            }
            "library_entry_id"
                if compact.strip_prefix("entry-").is_some_and(|hash| {
                    hash.len() == 64 && hash.bytes().all(|byte| byte.is_ascii_hexdigit())
                }) =>
            {
                compact
            }
            "library_proposal_id"
                if compact.strip_prefix("annotation-").is_some_and(|hash| {
                    hash.len() == 64 && hash.bytes().all(|byte| byte.is_ascii_hexdigit())
                }) =>
            {
                compact
            }
            "oligo_sequence" | "annealing_sequence" | "5prime_tail"
                if !compact.is_empty() && Sequence::dna(&compact).is_ok() =>
            {
                compact
            }
            _ => raw.clone(),
        }
    });
    (key.into(), normalized)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ImportReport {
    pub molecules: Vec<Molecule>,
    pub warnings: Vec<String>,
}

/// Content-based file import. Trace imports retain exact original attachments in native metadata.
pub fn import_bytes(bytes: &[u8]) -> Result<ImportReport, Error> {
    if bytes.len() > MAX_IMPORT_BYTES {
        return Err(Error::TooLarge);
    }
    if crate::legacy_formats::is_collagene_xml(bytes) {
        return crate::legacy_formats::import_collagene_xml(bytes);
    }
    if crate::legacy_formats::is_xdna(bytes) {
        return crate::legacy_formats::import_xdna(bytes, "Imported Strider document");
    }
    if bytes.len() >= 13 && &bytes[5..13] == b"SnapGene" {
        return crate::snapgene::import_snapgene(bytes);
    }
    if bytes.starts_with(b"ABIF") || bytes.starts_with(b".scf") {
        let trace = crate::traces::import_trace(bytes)?;
        return Ok(ImportReport {
            molecules: vec![trace.attached_molecule()?],
            warnings: trace.warnings,
        });
    }
    if bytes.starts_with(b"CG2\0") {
        return Ok(ImportReport {
            molecules: vec![crate::protocol::load_document(bytes)?],
            warnings: vec![],
        });
    }
    let text = std::str::from_utf8(bytes)
        .map_err(|_| Error::Unsupported("unrecognized binary sequence format".into()))?;
    import_text(text)
}

const PROVENANCE_KEYS: [&str; 12] = [
    "source:provider",
    "source:canonical_url",
    "source:identifier",
    "source:accession_version",
    "source:retrieved_at",
    "source:content_blake3",
    "source:format",
    "source:sequence_id",
    "source:addgene_origin",
    "source:addgene_completeness",
    "source:addgene_catalog_blake3",
    "source:addgene_sequence_description",
];
/// Safe, bounded source data for annotated text exports; credentials never belong in these fields.
pub fn provenance_comments(molecule: &Molecule) -> Vec<String> {
    let mut comments: Vec<String> = PROVENANCE_KEYS
        .iter()
        .filter_map(|key| {
            molecule.metadata.get(*key).map(|value| {
                format!(
                    "Collagene provenance {key}={}",
                    value
                        .chars()
                        .filter(|c| !c.is_control())
                        .take(8192)
                        .collect::<String>()
                )
            })
        })
        .collect();
    // These are portable producer claims, never imported scientific authority.
    let recipe = molecule
        .metadata
        .get(crate::restriction_inspection::RECIPE_KEY)
        .or_else(|| {
            molecule
                .metadata
                .get("source:restriction_inspection:recipe")
        });
    let ends = molecule.metadata.get("duplex_model").or_else(|| {
        molecule
            .metadata
            .get("source:restriction_inspection:duplex_summary")
    });
    if recipe.is_some() || ends.is_some() {
        comments.push(DIGEST_PRODUCER_MARKER.into());
        comments.push("Complete enzyme geometry, options and end summaries are unverified producer data. Sequence text cannot restore strand/end chemistry or a passing scientific check; retain native source-event ancestry and recompute.".into());
        if let Some(value) = recipe {
            digest_claim_comments(&mut comments, "CG-DIGEST-RECIPE-v1", value);
        }
        if let Some(value) = ends {
            digest_claim_comments(&mut comments, "CG-DIGEST-ENDS-v1", value);
        }
    }
    // Limits cover the full supported codon profile, including 10,000 residues
    // and its 32 KiB options; larger manually inserted values disclose truncation.
    let codon_keys = [
        ("codon_design:engine", 128),
        ("codon_design:origin", 256),
        ("codon_design:source_revision", 64),
        ("codon_design:options_revision", 64),
        ("codon_design:required_protein", 10_000),
        ("codon_design:options", 32 * 1024),
    ];
    if codon_keys
        .iter()
        .any(|(key, _)| molecule.metadata.contains_key(*key))
    {
        comments.push("Collagene codon design: unverified producer evidence; retain the source and report and recompute before use.".into());
        for (key, limit) in codon_keys {
            if let Some(value) = molecule.metadata.get(key) {
                let mut chars = value.chars().filter(|c| !c.is_control());
                let text: String = chars.by_ref().take(limit).collect();
                let truncated = if chars.next().is_some() {
                    " [truncated; retain native evidence]"
                } else {
                    ""
                };
                comments.push(format!("Collagene provenance {key}={text}{truncated}"));
            }
        }
    }
    let contig_keys = [
        ("contig:engine", 128),
        ("contig:batch_revision", 64),
        ("contig:options_revision", 64),
        ("contig:evidence_hash", 64),
        ("contig:scientific_product_revision", 64),
        ("contig:id", 32),
        ("contig:complete", 8),
        ("contig:layout_unique", 8),
        ("contig:options", 32 * 1024),
        ("contig:reads", 32 * 1024),
    ];
    if contig_keys
        .iter()
        .any(|(key, _)| molecule.metadata.contains_key(*key))
    {
        comments.push("Collagene read assembly: unverified producer evidence; completeness refers to the declared computational model, not biological identity. Retain source reads and native column evidence.".into());
        for (key, limit) in contig_keys {
            if let Some(value) = molecule.metadata.get(key) {
                let mut chars = value.chars().filter(|c| !c.is_control());
                let text: String = chars.by_ref().take(limit).collect();
                let suffix = if chars.next().is_some() {
                    " [truncated; retain native evidence]"
                } else {
                    ""
                };
                comments.push(format!("Collagene provenance {key}={text}{suffix}"));
            }
        }
    }
    if let Some(value) = molecule.metadata.get("tracked_modifications") {
        comments.push("Collagene positional modification summary: unverified producer evidence, not measured chemical state. Reconstruct the original plan and source assumptions before use; sequence text does not preserve strand ledgers or current verification.".into());
        let mut chars = value.chars().filter(|c| !c.is_control());
        let text: String = chars.by_ref().take(16 * 1024).collect();
        let suffix = if chars.next().is_some() {
            " [truncated; retain native graph]"
        } else {
            ""
        };
        comments.push(format!(
            "Collagene provenance tracked_modifications={text}{suffix}"
        ));
    }
    let edit_keys = [
        ("edit_template:engine", 128),
        ("edit_template:role", 32),
        ("edit_template:source_revision", 64),
        ("edit_template:target_revision", 64),
        ("edit_template:options_revision", 64),
        ("edit_template:design_revision", 64),
        ("edit_template:options", 32 * 1024),
        ("edit_template:guide", 16 * 1024),
        ("crispr:source_document", 128),
        ("crispr:source_event", 128),
        ("crispr:target_document", 128),
        ("crispr:target_event", 128),
        ("crispr:target_revision", 64),
        ("crispr:scientific_product_revision", 64),
    ];
    if edit_keys
        .iter()
        .any(|(key, _)| molecule.metadata.contains_key(*key))
    {
        comments.push("Collagene CRISPR/edit-template design: unverified producer sequence evidence. Source and independent-target requirements must be checked afresh; no cellular editing, guide activity or experimental ancestry is established.".into());
        for (key, limit) in edit_keys {
            if let Some(value) = molecule.metadata.get(key) {
                let mut chars = value.chars().filter(|c| !c.is_control());
                let text: String = chars.by_ref().take(limit).collect();
                let suffix = if chars.next().is_some() {
                    " [truncated; retain native source/target and full JSON]"
                } else {
                    ""
                };
                comments.push(format!("Collagene provenance {key}={text}{suffix}"));
            }
        }
    }
    comments
}
#[derive(Debug, Clone, Copy)]
pub enum TextFormat {
    GenBank,
    Embl,
    Fasta,
}
#[derive(Debug, Clone)]
pub struct ExportReport {
    pub text: String,
    pub warnings: Vec<String>,
}
pub fn export_with_report(molecule: &Molecule, format: TextFormat) -> Result<ExportReport, Error> {
    let text = match format {
        TextFormat::GenBank => export_genbank(molecule)?,
        TextFormat::Embl => crate::embl::export_embl(molecule)?,
        TextFormat::Fasta => export_fasta(molecule),
    };
    let mut warnings = vec![];
    if matches!(format, TextFormat::Fasta) {
        warnings.push("FASTA exports sequence only; features, topology, source provenance and attachments require the native or annotated source format.".into());
    }
    if !matches!(format, TextFormat::Embl) && molecule.metadata.contains_key("embl.original_record")
    {
        warnings.push("Original EMBL headers/evidence are retained in native/EMBL exports; this export contains only mapped semantics.".into());
    }
    if !matches!(format, TextFormat::GenBank) && molecule.original_genbank.is_some() {
        warnings.push("Original GenBank header/evidence fields are not all mapped to this export; keep the native or original GenBank record.".into());
    }
    if molecule
        .metadata
        .contains_key(crate::snapgene::ORIGINAL_BASE64_KEY)
    {
        warnings.push("SnapGene primers/history/attachments and unmapped source details are not preserved by this text export. Keep .cg2 or the original .dna evidence.".into());
    }
    if molecule.metadata.contains_key("trace_original_base64") {
        warnings.push("Chromatogram signals, peaks and caller-accuracy evidence are omitted by this text export. Keep .cg2 or the original trace.".into());
    }
    if molecule
        .metadata
        .contains_key(crate::oligo_formats::FASTQ_ORIGINAL_KEY)
    {
        warnings.push("FASTQ per-base qualities, explicit encoding and original record evidence are omitted by this text export. Keep the native workspace or original FASTQ.".into());
    }
    if molecule
        .metadata
        .keys()
        .any(|key| key.starts_with("oligo:"))
    {
        warnings.push("Oligo source rows, mapping and quantity/unit evidence are not fully represented by this text export. Keep the native workspace and original oligo table.".into());
    }
    if molecule
        .metadata
        .keys()
        .any(|key| key.starts_with("feature_library:"))
    {
        warnings.push("Feature-library identifiers are provenance references. Complete library templates/source evidence and reconstructed search reports are not represented by this text export; retain the workspace/library JSON and search report. Inferred annotations still require review.".into());
    }
    if molecule
        .metadata
        .keys()
        .any(|key| key.starts_with("annotation_import:"))
    {
        warnings.push("Original GFF3/BED annotation tables, explicit reference mappings and import applications require native storage or the original source. Rich source-row qualifiers may not survive INSDC physical line wrapping exactly; this text export is not a lossless annotation-dataset interchange.".into());
    }
    if molecule
        .metadata
        .contains_key(crate::legacy_formats::ORIGINAL_BASE64_KEY)
    {
        warnings.push("Original legacy XML/XDNA source fields and primer/aligned-strand evidence are not fully represented by this text export; retain the native workspace or original legacy file.".into());
    }
    if molecule
        .metadata
        .keys()
        .any(|key| key.starts_with("codon_design:"))
    {
        warnings.push("Codon-design comments are bounded, unverified producer evidence. Native storage and the JSON report retain the source/options bindings; text import does not restore a passing check or computational event parents.".into());
    }
    if molecule
        .metadata
        .keys()
        .any(|key| key.starts_with("contig:"))
    {
        warnings.push("Read-assembly comments are unverified producer evidence. Per-call quality, gaps, overlap alternatives, source events and current replay authority require the native workspace and full assembly report; this text export is not a lossless evidence interchange.".into());
    }
    if molecule.metadata.contains_key("tracked_modifications") {
        warnings.push("Positional modification comments contain only bounded, unverified producer summaries. Exact strand state and source-generation assumptions require fresh reconstruction from the native plan/source graph; this sequence export cannot restore chemical or verification authority.".into());
    }
    if molecule
        .metadata
        .keys()
        .any(|key| key.starts_with("edit_template:") || key.starts_with("crispr:"))
    {
        warnings.push("CRISPR/edit-template text comments are bounded, unverified producer summaries. Exact source/target event requirements, full guide/options/report JSON and current checks require native storage and fresh computation; FASTA and sequence text are not lossless design or verification interchange.".into());
    }
    if molecule.metadata.contains_key("duplex_model") {
        warnings.push("This is a reference-strand sequence export. Duplex end states/nicks are report evidence and must be reconstructed from their originating plan.".into());
    }
    if molecule
        .metadata
        .contains_key(crate::restriction_inspection::RECIPE_KEY)
        || molecule
            .metadata
            .contains_key("source:restriction_inspection:recipe")
    {
        warnings.push("Digest recipe/end comments are bounded, unverified producer data. This reference-strand export cannot restore physical end states or source-event ancestry; retain the native workspace and recompute from its exact source. FASTA is sequence-only.".into());
    }
    Ok(ExportReport { text, warnings })
}

pub fn import_text(text: &str) -> Result<ImportReport, Error> {
    if crate::legacy_formats::is_collagene_xml(text.as_bytes()) {
        return crate::legacy_formats::import_collagene_xml(text.as_bytes());
    }

    if text.len() > MAX_IMPORT_BYTES {
        return Err(Error::TooLarge);
    }
    let text = text.trim_start_matches('\u{feff}').trim_start();
    if text.starts_with("LOCUS") {
        import_genbank(text)
    } else if text.starts_with("ID ") {
        crate::embl::import_embl(text)
    } else if text.starts_with('>') {
        import_fasta(text)
    } else {
        let mol = Molecule::dna("Untitled sequence", text, Topology::Linear)?;
        Ok(ImportReport {
            molecules: vec![mol],
            warnings: vec![
                "Plain sequence has no topology or annotations; imported as linear DNA.".into(),
            ],
        })
    }
}

pub fn import_fasta(text: &str) -> Result<ImportReport, Error> {
    if text.len() > MAX_IMPORT_BYTES {
        return Err(Error::TooLarge);
    }
    let mut molecules = vec![];
    let mut name: Option<String> = None;
    let mut bases = String::new();
    let mut alphabet = Alphabet::Dna;
    let mut warnings = vec!["FASTA does not carry feature annotations or explicit topology; records are imported as linear.".into()];
    let finish = |name: String, bases: &str, alphabet| -> Result<Molecule, Error> {
        let mut mol = Molecule::dna(name, "N", Topology::Linear)?;
        mol.sequence = Sequence::new(alphabet, bases)?;
        Ok(mol)
    };
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if let Some(header) = line.strip_prefix('>') {
            if let Some(name) = name.take() {
                molecules.push(finish(name, &bases, alphabet)?);
            }
            if molecules.len() >= 1000 {
                return Err(Error::TooLarge);
            }
            if header.trim().is_empty() {
                return Err(Error::Format("FASTA header is empty".into()));
            }
            name = Some(header.trim().into());
            bases.clear();
            // Explicit annotation avoids guessing DNA versus protein for ambiguous short records.
            alphabet = if header.contains("[molecule=protein]") {
                Alphabet::Protein
            } else if header.contains("[molecule=RNA]") {
                Alphabet::Rna
            } else {
                Alphabet::Dna
            };
        } else {
            if name.is_none() {
                return Err(Error::Format(
                    "FASTA sequence appears before its header".into(),
                ));
            }
            if bases.len() + line.len() > MAX_SEQUENCE_LEN {
                return Err(Error::TooLarge);
            }
            bases.push_str(line);
        }
    }
    if let Some(name) = name {
        molecules.push(finish(name, &bases, alphabet)?);
    }
    if molecules.is_empty() {
        return Err(Error::Format("no FASTA records".into()));
    }
    if molecules
        .iter()
        .any(|m| m.sequence.alphabet() != Alphabet::Dna)
    {
        warnings.push(
            "Non-DNA alphabet was taken from an explicit [molecule=...] FASTA header.".into(),
        );
    }
    Ok(ImportReport {
        molecules,
        warnings,
    })
}

fn check_location_text(text: &str) -> Result<(), Error> {
    if text.is_empty()
        || text.len() > 65_536
        || !text.is_ascii()
        || text.bytes().any(|b| b.is_ascii_control())
    {
        return Err(Error::Format(
            "invalid or oversized location expression".into(),
        ));
    }
    let mut depth = 0usize;
    for b in text.bytes() {
        match b {
            b'(' => {
                depth += 1;
                if depth >= 32 {
                    return Err(Error::Format("location nesting exceeds 32 levels".into()));
                }
            }
            b')' => {
                depth = depth
                    .checked_sub(1)
                    .ok_or_else(|| Error::Format("unbalanced location parentheses".into()))?;
            }
            _ => {}
        }
    }
    if depth != 0 {
        return Err(Error::Format("unbalanced location parentheses".into()));
    }
    Ok(())
}
pub fn parse_location(text: &str, length: usize) -> Result<Location, Error> {
    let text = text.trim();
    check_location_text(text)?;
    let raw = gb::Location::from_gb_format(text)
        .map_err(|e| Error::Format(format!("invalid location: {e}")))?;
    let location = from_location(&raw, 0)?;
    location.validate(length)?;
    Ok(location)
}
pub fn location_text(location: &Location) -> Result<String, Error> {
    location.validate(MAX_SEQUENCE_LEN)?;
    if let Location::Unresolved(raw) = location {
        check_location_text(raw)?;
        return Ok(raw.clone());
    }
    Ok(to_location(location)?.to_string())
}

fn from_location(l: &gb::Location, depth: usize) -> Result<Location, Error> {
    if depth >= 32 {
        return Err(Error::Format("location nesting exceeds 32 levels".into()));
    }
    Ok(match l {
        gb::Location::Range((a, before), (b, after)) => Location::Range {
            start: usize::try_from(*a)
                .map_err(|_| Error::Coordinates("negative GenBank start".into()))?,
            end: usize::try_from(*b)
                .map_err(|_| Error::Coordinates("negative GenBank end".into()))?,
            before: before.0,
            after: after.0,
        },
        gb::Location::Complement(inner) => {
            Location::Complement(Box::new(from_location(inner, depth + 1)?))
        }
        gb::Location::Join(parts) => Location::Join(
            parts
                .iter()
                .map(|p| from_location(p, depth + 1))
                .collect::<Result<_, _>>()?,
        ),
        gb::Location::Order(parts) => Location::Order(
            parts
                .iter()
                .map(|p| from_location(p, depth + 1))
                .collect::<Result<_, _>>()?,
        ),
        _ => Location::Unresolved(l.to_string()),
    })
}

fn to_location(l: &Location) -> Result<gb::Location, Error> {
    Ok(match l {
        Location::Range {
            start,
            end,
            before,
            after,
        } => gb::Location::Range(
            ((*start) as i64, Before(*before)),
            ((*end) as i64, After(*after)),
        ),
        Location::Complement(inner) => gb::Location::Complement(Box::new(to_location(inner)?)),
        Location::Join(parts) => {
            gb::Location::Join(parts.iter().map(to_location).collect::<Result<_, _>>()?)
        }
        Location::Order(parts) => {
            gb::Location::Order(parts.iter().map(to_location).collect::<Result<_, _>>()?)
        }
        Location::Unresolved(raw) => {
            check_location_text(raw)?;
            gb::Location::from_gb_format(raw)
                .map_err(|e| Error::Format(format!("invalid retained location: {e}")))?
        }
    })
}

/// Bound recursive location syntax before invoking the general GenBank parser.
fn preflight_genbank(text: &str) -> Result<(), Error> {
    let mut in_features = false;
    let mut in_location = false;
    let mut depth = 0usize;
    let mut features = 0usize;
    for line in text.lines() {
        if line.starts_with("FEATURES") {
            in_features = true;
            continue;
        }
        if line.starts_with("ORIGIN") || line.starts_with("CONTIG") || line.starts_with("//") {
            in_features = false;
        }
        if !in_features {
            continue;
        }
        let Some(key) = line.get(5..21) else {
            continue;
        };
        let Some(value) = line.get(21..) else {
            continue;
        };
        if !key.trim().is_empty() {
            features += 1;
            if features > 100_000 {
                return Err(Error::TooLarge);
            }
            in_location = true;
            depth = 0;
        }
        if value.trim_start().starts_with('/') {
            in_location = false;
        }
        if in_location {
            for c in value.chars() {
                if c == '(' {
                    depth += 1;
                    if depth >= 32 {
                        return Err(Error::Format(
                            "GenBank location nesting exceeds 32 levels".into(),
                        ));
                    }
                } else if c == ')' {
                    depth = depth.saturating_sub(1);
                }
            }
        }
    }
    Ok(())
}

pub fn import_genbank(text: &str) -> Result<ImportReport, Error> {
    if text.len() > MAX_IMPORT_BYTES {
        return Err(Error::TooLarge);
    }
    preflight_genbank(text)?;
    // Retain each record once. Copying the full multi-record input into every molecule
    // would multiply memory use and make native exports unnecessarily enormous.
    let mut raw_records = vec![];
    let mut start = 0usize;
    let mut offset = 0usize;
    for line in text.split_inclusive('\n') {
        offset += line.len();
        if line.trim_end() == "//" && line.starts_with("//") {
            raw_records.push(&text[start..offset]);
            start = offset;
            if raw_records.len() > 1000 {
                return Err(Error::TooLarge);
            }
        }
    }
    if start < text.len() && !text[start..].trim().is_empty() {
        raw_records.push(&text[start..]);
    }
    let mut molecules = vec![];
    let mut warnings = vec![];
    for (index, record) in gb_io::reader::SeqReader::new(text.as_bytes()).enumerate() {
        if index >= 1000 {
            return Err(Error::TooLarge);
        }
        let record = record.map_err(|e| Error::Format(format!("GenBank: {e}")))?;
        if record.seq.len() > MAX_SEQUENCE_LEN {
            return Err(Error::TooLarge);
        }
        if record.seq.is_empty() {
            return Err(Error::Unsupported(
                "GenBank CONTIG-only records require reference resolution".into(),
            ));
        }
        if record.len.is_some_and(|len| len != record.seq.len()) {
            return Err(Error::Format(
                "GenBank LOCUS length differs from sequence length".into(),
            ));
        }
        let molecule_type = record
            .molecule_type
            .as_deref()
            .unwrap_or("DNA")
            .to_ascii_uppercase();
        let alphabet = if molecule_type.contains("RNA") {
            Alphabet::Rna
        } else if molecule_type == "AA" || molecule_type.contains("PROTEIN") {
            Alphabet::Protein
        } else {
            Alphabet::Dna
        };
        let bases = std::str::from_utf8(&record.seq)
            .map_err(|_| Error::Format("non-ASCII sequence".into()))?;
        let mut mol = Molecule::dna(
            record
                .name
                .clone()
                .unwrap_or_else(|| format!("Sequence {}", index + 1)),
            "N",
            Topology::Linear,
        )?;
        mol.sequence = Sequence::new(alphabet, bases)?;
        if let Some(kind) = &record.molecule_type {
            mol.metadata
                .insert("genbank.molecule_type".into(), kind.clone());
        }

        mol.topology = if record.is_circular() {
            Topology::Circular
        } else {
            Topology::Linear
        };
        mol.original_genbank = Some(
            raw_records
                .get(index)
                .ok_or_else(|| Error::Format("cannot associate original GenBank record".into()))?
                .to_string(),
        );
        if retain_digest_producer_comments(&mut mol, raw_records[index]) {
            warnings.push("Sequence-only digest producer comments were retained as unverified source data; strand/end chemistry was not restored.".into());
        }
        mol.metadata
            .insert("source_record_index".into(), "0".into());
        if let Some(def) = &record.definition {
            mol.metadata.insert("definition".into(), def.clone());
        }
        for f in &record.features {
            let location = from_location(&f.location, 0)?;
            if has_unresolved(&location) {
                warnings.push(format!(
                    "{}: {} location retained but not interpreted: {}",
                    mol.name, f.kind, f.location
                ));
            }
            mol.features.push(Feature {
                kind: f.kind.to_string(),
                location,
                qualifiers: f
                    .qualifiers
                    .iter()
                    .map(|(k, v)| import_qualifier(k, v))
                    .collect(),
                needs_review: f.qualifiers.iter().any(|(key, value)| {
                    key == "collagene_review_required" && value.as_deref() == Some("true")
                }),
            });
        }
        mol.validate()?;
        molecules.push(mol);
    }
    if molecules.is_empty() {
        return Err(Error::Format("no GenBank records".into()));
    }
    warnings.push("Original GenBank text is retained. Headers/qualifiers supported by the parser are exported; unrecognized syntax may require the original attachment.".into());
    Ok(ImportReport {
        molecules,
        warnings,
    })
}

fn has_unresolved(l: &Location) -> bool {
    match l {
        Location::Unresolved(_) => true,
        Location::Complement(v) => has_unresolved(v),
        Location::Join(v) | Location::Order(v) => v.iter().any(has_unresolved),
        _ => false,
    }
}

pub fn export_genbank(mol: &Molecule) -> Result<String, Error> {
    mol.validate()?;
    if mol.sequence.alphabet() == Alphabet::Protein {
        return Err(Error::Unsupported(
            "protein GenPept export is not implemented; use FASTA or the native format".into(),
        ));
    }
    let mut record = if let Some(original) = &mol.original_genbank {
        let index = mol
            .metadata
            .get("source_record_index")
            .and_then(|i| i.parse::<usize>().ok())
            .unwrap_or(0);
        gb_io::reader::SeqReader::new(original.as_bytes())
            .nth(index)
            .ok_or_else(|| Error::Format("original record index is missing".into()))?
            .map_err(|e| Error::Format(e.to_string()))?
    } else {
        gb::Seq::empty()
    };
    if !record.references.is_empty()
        && !record
            .seq
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .eq(mol.sequence.as_str().bytes())
    {
        return Err(Error::Unsupported("GenBank reference ranges cannot yet be remapped after sequence changes; use the native format to retain the original evidence".into()));
    }
    record.name = Some(mol.name.replace(char::is_whitespace, "_"));
    record.seq = mol.sequence.as_str().as_bytes().to_vec();
    record.len = Some(record.seq.len());
    record.contig = None;
    record.topology = if mol.topology == Topology::Circular {
        gb::Topology::Circular
    } else {
        gb::Topology::Linear
    };
    let saved_type = mol.metadata.get("genbank.molecule_type").filter(|value| {
        value.len() <= 32
            && value
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || b == b'-')
    });
    let saved_type = saved_type.filter(|value| match mol.sequence.alphabet() {
        Alphabet::Dna => value.to_ascii_uppercase().ends_with("DNA"),
        Alphabet::Rna => value.to_ascii_uppercase().ends_with("RNA"),
        Alphabet::Protein => false,
    });
    record.molecule_type = Some(saved_type.cloned().unwrap_or_else(|| {
        match mol.sequence.alphabet() {
            Alphabet::Dna => "DNA",
            Alphabet::Rna => "RNA",
            Alphabet::Protein => "AA",
        }
        .into()
    }));
    record.definition = mol.metadata.get("definition").cloned();
    for comment in provenance_comments(mol) {
        if !record.comments.contains(&comment) {
            record.comments.push(comment);
        }
    }
    record.features = mol
        .features
        .iter()
        .map(|f| {
            let mut qualifiers: Vec<_> = f
                .qualifiers
                .iter()
                .filter(|(key, _)| key != "collagene_review_required")
                .map(|(k, v)| (Cow::Owned(k.clone()), v.clone()))
                .collect();
            if f.needs_review {
                qualifiers.push((
                    Cow::Borrowed("collagene_review_required"),
                    Some("true".into()),
                ));
                if !qualifiers.iter().any(|(key, value)| {
                    key == "note"
                        && value.as_deref() == Some("Collagene: annotation requires review")
                }) {
                    qualifiers.push((
                        Cow::Borrowed("note"),
                        Some("Collagene: annotation requires review".into()),
                    ));
                }
            }
            Ok(gb::Feature {
                kind: Cow::Owned(f.kind.clone()),
                location: to_location(&f.location)?,
                qualifiers,
            })
        })
        .collect::<Result<_, Error>>()?;
    let mut bytes = vec![];
    record
        .write(&mut bytes)
        .map_err(|e| Error::Format(e.to_string()))?;
    String::from_utf8(bytes).map_err(|e| Error::Format(e.to_string()))
}

pub fn export_fasta(mol: &Molecule) -> String {
    let alphabet = match mol.sequence.alphabet() {
        Alphabet::Rna if !mol.name.contains("[molecule=RNA]") => " [molecule=RNA]",
        Alphabet::Protein if !mol.name.contains("[molecule=protein]") => " [molecule=protein]",
        _ => "",
    };
    let mut out = format!(">{}{alphabet}\n", mol.name.replace(['\n', '\r'], " "));
    for chunk in mol.sequence.as_str().as_bytes().chunks(80) {
        out.push_str(std::str::from_utf8(chunk).expect("ASCII"));
        out.push('\n');
    }
    out
}
