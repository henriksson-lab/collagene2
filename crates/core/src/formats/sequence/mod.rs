//! GenBank syntax uses the standalone gb-io parser; molecular semantics remain ours.
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_model::*;
use gb_io::seq::{self as gb, After, Before};
use serde::{Deserialize, Serialize};
use std::borrow::Cow;

pub const MAX_IMPORT_BYTES: usize = 32 * 1024 * 1024;
const DIGEST_PRODUCER_MARKER: &str =
    "Collagene physical product v1: unverified sequence-only evidence";
const REFERENCE_PRODUCER_MARKER: &str =
    "Collagene reference-only v1: unverified sequence-only evidence";
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
    let reference_only = raw
        .lines()
        .filter_map(comment)
        .any(|s| s.trim_end() == REFERENCE_PRODUCER_MARKER);
    if !reference_only
        && !raw
            .lines()
            .filter_map(comment)
            .any(|s| s.trim_end() == DIGEST_PRODUCER_MARKER)
    {
        return false;
    }
    m.metadata
        .insert("source:legacy:chemistry_unmapped".into(), "true".into());
    if reference_only {
        m.metadata
            .insert("source:recombination:reference_only".into(), "true".into());
    }
    for (prefix, key) in [
        (
            "CG-RECOMB-REF-v1 ",
            "source:recombination:unverified_recipe",
        ),
        (
            "CG-DIGEST-RECIPE-v1 ",
            "source:restriction_inspection:recipe",
        ),
        (
            "CG-DIGEST-ENDS-v1 ",
            "source:restriction_inspection:duplex_summary",
        ),
        (
            "CG-CATALOG-DIGEST-v1 ",
            "source:restriction_catalog:selection",
        ),
        ("CG-CATALOG-ENGINE-v1 ", "source:restriction_catalog:engine"),
        (
            "CG-CATALOG-REVISION-v1 ",
            "source:restriction_catalog:batch_revision",
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
    if crate::formats::original::is_collagene_xml(bytes) {
        return crate::formats::original::import_collagene_xml(bytes);
    }
    if crate::formats::original::is_xdna(bytes) {
        return crate::formats::original::import_xdna(bytes, "Imported Strider document");
    }
    if bytes.len() >= 13 && &bytes[5..13] == b"SnapGene" {
        return crate::formats::snapgene::import_snapgene(bytes);
    }
    if bytes.starts_with(b"ABIF") || bytes.starts_with(b".scf") {
        let trace = crate::review::reads::traces::import_trace(bytes)?;
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
        .get(crate::cloning::restriction::inspection::RECIPE_KEY)
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
    let catalog = molecule
        .metadata
        .get(crate::cloning::restriction::catalog::RECIPE_KEY)
        .or_else(|| {
            molecule
                .metadata
                .get("source:restriction_catalog:selection")
        });
    let catalog_engine = molecule
        .metadata
        .get(crate::cloning::restriction::catalog::ENGINE_KEY)
        .or_else(|| molecule.metadata.get("source:restriction_catalog:engine"));
    let catalog_revision = molecule
        .metadata
        .get(crate::cloning::restriction::catalog::PRODUCT_REVISION_KEY)
        .or_else(|| {
            molecule
                .metadata
                .get("source:restriction_catalog:batch_revision")
        });
    let reference_recipe = molecule
        .metadata
        .get(crate::cloning::recombination::RECIPE_KEY)
        .or_else(|| {
            molecule
                .metadata
                .get("source:recombination:unverified_recipe")
        });
    let reference_only = molecule
        .metadata
        .get(crate::cloning::recombination::REFERENCE_ONLY_KEY)
        .is_some_and(|s| s == "true")
        || molecule
            .metadata
            .get("source:recombination:reference_only")
            .is_some_and(|s| s == "true")
        || reference_recipe.is_some();
    let unmapped = molecule
        .metadata
        .get("source:legacy:chemistry_unmapped")
        .is_some_and(|s| s == "true");
    if reference_only {
        comments.push(REFERENCE_PRODUCER_MARKER.into());
        comments.push("Reference recombination recipes are unverified producer claims. These bases do not establish physical strand states, methylation, recombinase activity or cellular repair. Retain native source/profile history and reconstruct afresh.".into());
        if let Some(value) = reference_recipe {
            digest_claim_comments(&mut comments, "CG-RECOMB-REF-v1", value);
        }
    }
    if recipe.is_some()
        || ends.is_some()
        || catalog.is_some()
        || catalog_engine.is_some()
        || catalog_revision.is_some()
        || (unmapped && !reference_only)
    {
        comments.push(DIGEST_PRODUCER_MARKER.into());
        comments.push("Complete enzyme geometry, options and end summaries are unverified producer data. Sequence text cannot restore strand/end chemistry or a passing scientific check; retain native source-event ancestry and recompute.".into());
        if let Some(value) = recipe {
            digest_claim_comments(&mut comments, "CG-DIGEST-RECIPE-v1", value);
        }
        if let Some(value) = ends {
            digest_claim_comments(&mut comments, "CG-DIGEST-ENDS-v1", value);
        }
        if let Some(value) = catalog {
            comments.push("Catalog edition, content hash, ordered enzyme identities and geometry are unverified producer claims. Reimport does not authenticate supplier evidence or restore current catalog, cleavage or verification authority.".into());
            digest_claim_comments(&mut comments, "CG-CATALOG-DIGEST-v1", value);
        }
        if let Some(value) = catalog_engine {
            digest_claim_comments(&mut comments, "CG-CATALOG-ENGINE-v1", value);
        }
        if let Some(value) = catalog_revision {
            digest_claim_comments(&mut comments, "CG-CATALOG-REVISION-v1", value);
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
    GenPept,
    UniProt,
}
#[derive(Debug, Clone)]
pub struct ExportReport {
    pub text: String,
    pub warnings: Vec<String>,
}
pub fn export_with_report(molecule: &Molecule, format: TextFormat) -> Result<ExportReport, Error> {
    let mut warnings = vec![];
    let text = match format {
        TextFormat::GenPept => {
            let report = crate::formats::protein::export(
                std::slice::from_ref(molecule),
                crate::formats::protein::Profile::GenPept,
            )?;
            warnings.extend(report.warnings);
            report.text
        }
        TextFormat::UniProt => {
            let report = crate::formats::protein::export(
                std::slice::from_ref(molecule),
                crate::formats::protein::Profile::UniProt,
            )?;
            warnings.extend(report.warnings);
            report.text
        }
        TextFormat::GenBank => {
            if molecule.sequence.alphabet() == Alphabet::Rna {
                warnings.push("RNA is stored internally as U; GenBank ORIGIN uses the INSDC T encoding. The native document and retained original remain unchanged.".into());
            }
            export_genbank(molecule)?
        }
        TextFormat::Embl => crate::formats::sequence::embl::export_embl(molecule)?,
        TextFormat::Fasta => export_fasta(molecule),
    };
    if matches!(format, TextFormat::Fasta) {
        warnings.push("FASTA exports sequence only; features, topology, source provenance and attachments require the native or annotated source format.".into());
        if fasta_export_header(molecule).1 {
            warnings.push("Conflicting alphabet markers in the literal document name were escaped in the FASTA header; the emitted canonical marker describes the actual sequence alphabet. The native name remains unchanged.".into());
        }
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
        .contains_key(crate::formats::snapgene::ORIGINAL_BASE64_KEY)
    {
        warnings.push("SnapGene primers/history/attachments and unmapped source details are not preserved by this text export. Keep .cg2 or the original .dna evidence.".into());
    }
    if molecule.metadata.contains_key("trace_original_base64") {
        warnings.push("Chromatogram signals, peaks and caller-accuracy evidence are omitted by this text export. Keep .cg2 or the original trace.".into());
    }
    if molecule
        .metadata
        .contains_key(crate::formats::oligos::FASTQ_ORIGINAL_KEY)
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
        .contains_key(crate::formats::original::ORIGINAL_BASE64_KEY)
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
        .contains_key(crate::cloning::recombination::REFERENCE_ONLY_KEY)
        || molecule
            .metadata
            .contains_key("source:recombination:reference_only")
        || molecule
            .metadata
            .contains_key(crate::cloning::recombination::RECIPE_KEY)
        || molecule
            .metadata
            .contains_key("source:recombination:unverified_recipe")
        || molecule
            .metadata
            .get("source:legacy:chemistry_unmapped")
            .is_some_and(|s| s == "true")
    {
        warnings.push("Reference-only or unmapped producer evidence cannot restore physical chemistry on text reimport. Recombination recipes are bounded unverified claims; retain native source/profile history and reconstruct afresh. FASTA is sequence-only.".into());
    }
    if molecule
        .metadata
        .contains_key(crate::cloning::restriction::inspection::RECIPE_KEY)
        || molecule
            .metadata
            .contains_key("source:restriction_inspection:recipe")
    {
        warnings.push("Digest recipe/end comments are bounded, unverified producer data. This reference-strand export cannot restore physical end states or source-event ancestry; retain the native workspace and recompute from its exact source. FASTA is sequence-only.".into());
    }
    if molecule
        .metadata
        .contains_key(crate::cloning::restriction::catalog::RECIPE_KEY)
        || molecule
            .metadata
            .contains_key("source:restriction_catalog:selection")
    {
        warnings.push("Catalog-pinned digest comments retain bounded unverified edition/hash/identity/geometry producer claims only. Text import does not restore supplier facts, physical chemistry or fresh catalog matching; retain native source history and recompute. FASTA is sequence-only.".into());
    }
    Ok(ExportReport { text, warnings })
}

pub fn import_text(text: &str) -> Result<ImportReport, Error> {
    if crate::formats::original::is_collagene_xml(text.as_bytes()) {
        return crate::formats::original::import_collagene_xml(text.as_bytes());
    }

    if text.len() > MAX_IMPORT_BYTES {
        return Err(Error::TooLarge);
    }
    if crate::formats::sequence::gcg::is_gcg(text.as_bytes()) {
        return Err(Error::Unsupported("GCG single-sequence import requires an explicit DNA, RNA or protein alphabet; use the GCG import profile.".into()));
    }
    if let Some(profile) = crate::formats::protein::detect(text.as_bytes()) {
        return crate::formats::protein::import_bytes(text.as_bytes(), profile);
    }
    let text = text.trim_start_matches('\u{feff}').trim_start();
    if text.starts_with("LOCUS") {
        import_genbank(text)
    } else if text.starts_with("ID ") {
        crate::formats::sequence::embl::import_embl(text)
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

const FASTA_ALPHABET_MARKERS: [(&str, Alphabet); 3] = [
    ("[molecule=DNA]", Alphabet::Dna),
    ("[molecule=RNA]", Alphabet::Rna),
    ("[molecule=protein]", Alphabet::Protein),
];

fn fasta_header_alphabet(header: &str) -> Result<Option<Alphabet>, Error> {
    let mut declared = None;
    for (marker, alphabet) in FASTA_ALPHABET_MARKERS {
        if header.contains(marker) {
            if declared.is_some_and(|current| current != alphabet) {
                return Err(Error::Format(
                    "FASTA header contains conflicting explicit alphabet markers".into(),
                ));
            }
            declared = Some(alphabet);
        }
    }
    Ok(declared)
}

pub fn import_fasta(text: &str) -> Result<ImportReport, Error> {
    import_fasta_inner(text, None)
}

/// Explicit file/profile alphabet for unmarked FASTA, including protein strings
/// that also happen to be valid DNA ambiguity calls. No composition inference.
/// Existing recognized header declarations must agree with the selected profile.
pub fn import_fasta_with_alphabet(text: &str, alphabet: Alphabet) -> Result<ImportReport, Error> {
    import_fasta_inner(text, Some(alphabet))
}

fn import_fasta_inner(text: &str, explicit: Option<Alphabet>) -> Result<ImportReport, Error> {
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
            let declared = fasta_header_alphabet(header)?;
            if let (Some(selected), Some(declared)) = (explicit, declared) {
                if selected != declared {
                    return Err(Error::Format(
                        "FASTA header alphabet differs from the explicit import profile".into(),
                    ));
                }
            }
            alphabet = explicit.or(declared).unwrap_or(Alphabet::Dna);
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
    if let Some(alphabet) = explicit {
        warnings.push(format!(
            "FASTA alphabet was explicitly selected as {alphabet:?}; sequence composition and filename were not used to infer it."
        ));
    } else if molecules
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

// INSDC nucleotide flat files encode RNA uracil as T. Keep this conversion
// local to GenBank syntax: native RNA remains U, and original bytes stay exact.
fn genbank_alphabet(kind: Option<&str>) -> Alphabet {
    let kind = kind.unwrap_or("DNA").to_ascii_uppercase();
    if kind.contains("RNA") {
        Alphabet::Rna
    } else if kind == "AA" || kind.contains("PROTEIN") {
        Alphabet::Protein
    } else {
        Alphabet::Dna
    }
}
fn genbank_internal_base(base: u8, alphabet: Alphabet) -> u8 {
    let base = base.to_ascii_uppercase();
    if alphabet == Alphabet::Rna && base == b'T' {
        b'U'
    } else {
        base
    }
}
fn genbank_internal_sequence(bytes: &[u8], alphabet: Alphabet) -> Result<Sequence, Error> {
    let bases: Vec<u8> = bytes
        .iter()
        .map(|&b| genbank_internal_base(b, alphabet))
        .collect();
    let text = std::str::from_utf8(&bases)
        .map_err(|_| Error::Format("GenBank sequence is not ASCII".into()))?;
    Sequence::new(alphabet, text)
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
        let alphabet = genbank_alphabet(record.molecule_type.as_deref());
        let mut mol = Molecule::dna(
            record
                .name
                .clone()
                .unwrap_or_else(|| format!("Sequence {}", index + 1)),
            "N",
            Topology::Linear,
        )?;
        mol.sequence = genbank_internal_sequence(&record.seq, alphabet)?;
        if alphabet == Alphabet::Rna && record.seq.iter().any(|b| matches!(b, b't' | b'T')) {
            warnings.push(format!("{}: RNA-declared GenBank sequence normalized from T to U; original text retains the supplied encoding.", mol.name));
        }
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
            warnings.push("Sequence-only producer comments were retained as unverified source data; reference provenance and strand/end chemistry were not restored.".into());
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
    if mol
        .metadata
        .get(GENBANK_EDIT_RECIPE_KEY)
        .is_some_and(|raw| raw.len() > MAX_GENBANK_EDIT_RECIPE_BYTES)
    {
        return Err(Error::TooLarge);
    }
    mol.validate()?;
    if mol.sequence.alphabet() == Alphabet::Protein {
        return Err(Error::Unsupported(
            "protein documents require the explicit GenPept or UniProt export profile; FASTA and native exports are also available".into(),
        ));
    }
    let mut record = if let Some(original) = &mol.original_genbank {
        if original.len() > MAX_IMPORT_BYTES {
            return Err(Error::TooLarge);
        }
        // Bound the existing parser's LOCUS preallocation before rich parsing.
        if original
            .lines()
            .filter(|line| line.starts_with("LOCUS"))
            .any(|line| {
                line.split_whitespace()
                    .nth(2)
                    .and_then(|n| n.parse::<usize>().ok())
                    .is_some_and(|n| n > MAX_SEQUENCE_LEN)
            })
        {
            return Err(Error::TooLarge);
        }
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
    let original_alphabet = genbank_alphabet(record.molecule_type.as_deref());
    let source_differs = original_alphabet != mol.sequence.alphabet()
        || !record
            .seq
            .iter()
            .map(|&b| genbank_internal_base(b, original_alphabet))
            .eq(mol.sequence.as_str().bytes());
    if mol.metadata.contains_key(GENBANK_EDIT_RECIPE_KEY) {
        let index = mol
            .metadata
            .get("source_record_index")
            .map(|v| v.parse::<usize>())
            .transpose()
            .map_err(|_| Error::Unsupported("Invalid GenBank source record index".into()))?
            .unwrap_or(0);
        let original_sequence = genbank_internal_sequence(&record.seq, original_alphabet)?;
        let recipe = verified_genbank_edit_recipe(mol, &original_sequence, index)?
            .ok_or_else(|| Error::Unsupported("Missing GenBank edit recipe".into()))?;
        if record.references.len() > 1024 {
            return Err(Error::TooLarge);
        }
        let mut spans = 0usize;
        for reference in &record.references {
            if reference.description.len() > 4096 {
                return Err(Error::TooLarge);
            }
            spans = spans
                .checked_add(reference.description.bytes().filter(|&b| b == b';').count() + 1)
                .ok_or(Error::TooLarge)?;
        }
        if spans > 4096 || spans.saturating_mul(recipe.edits.len()) > 1_000_000 {
            return Err(Error::TooLarge);
        }
        for reference in &mut record.references {
            reference.description = remap_genbank_reference(&reference.description, &recipe)?;
            if !recipe.edits.is_empty() {
                let review = if mol.sequence.alphabet() == Alphabet::Rna {
                    "Collagene: source citation retained with projected coordinate scope; citation support for edited RNA is unverified. Original record remains immutable in native source history."
                } else {
                    "Collagene: source citation retained with projected coordinate scope; citation support for edited DNA is unverified. Original record remains immutable in native source history."
                };
                let remark = reference.remark.get_or_insert_with(String::new);
                if !remark.contains(review) {
                    if !remark.is_empty() {
                        remark.push(' ');
                    }
                    remark.push_str(review);
                }
            }
        }
    } else if !record.references.is_empty() && source_differs {
        return Err(Error::Unsupported("GenBank citation coordinates require a verified source-bound edit recipe; retain the original record or use the native format".into()));
    }
    if source_differs {
        for (kind, value) in [
            ("accession", record.accession.take()),
            ("version", record.version.take()),
        ] {
            if let Some(value) = value {
                if value.len() > 4096 {
                    return Err(Error::TooLarge);
                }
                record.comments.push(format!(
                    "Original source {kind} (ancestral provenance only): {value}"
                ));
            }
        }
    }
    record.name = Some(mol.name.replace(char::is_whitespace, "_"));
    record.seq = mol
        .sequence
        .as_str()
        .bytes()
        .map(|base| {
            if mol.sequence.alphabet() == Alphabet::Rna && base == b'U' {
                b'T'
            } else {
                base
            }
        })
        .collect();
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
    let (header, _) = fasta_export_header(mol);
    let mut out = format!(">{header}\n");
    for chunk in mol.sequence.as_str().as_bytes().chunks(80) {
        out.push_str(std::str::from_utf8(chunk).expect("ASCII"));
        out.push('\n');
    }
    out
}

fn fasta_export_header(mol: &Molecule) -> (String, bool) {
    let actual = mol.sequence.alphabet();
    let mut header = mol.name.replace(['\n', '\r'], " ");
    let mut escaped = false;
    for (marker, alphabet) in FASTA_ALPHABET_MARKERS {
        if alphabet != actual && header.contains(marker) {
            let label = format!(
                "{}; label]",
                marker.strip_suffix(']').expect("fixed marker")
            );
            header = header.replace(marker, &label);
            escaped = true;
        }
    }
    let marker = FASTA_ALPHABET_MARKERS
        .iter()
        .find(|(_, alphabet)| *alphabet == actual)
        .map(|(marker, _)| *marker)
        .expect("all native alphabets");
    // Keep established unmarked-DNA output, but disambiguate a name whose
    // conflicting semantic marker had to become a literal label.
    if (actual != Alphabet::Dna || escaped) && !header.contains(marker) {
        header.push(' ');
        header.push_str(marker);
    }
    (header, escaped)
}

pub mod embl;
pub mod export;
pub mod gcg;
