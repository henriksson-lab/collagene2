//! Official Addgene developer catalog adapter, independently implemented from its OpenAPI schema.
//! Signed download URLs and raw catalog responses are transient, never serialized as provenance.
use super::*;
use collagene_model::{Alphabet, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

pub(crate) const API_HOST: &str = "api.developers.addgene.org";
const MAX_CHOICES: usize = 1000;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum AddgeneSequenceOrigin {
    Depositor,
    Addgene,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum AddgeneCompleteness {
    Full,
    Partial,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct AddgeneSequenceChoice {
    pub sequence_id: u64,
    pub description: String,
    pub origin: AddgeneSequenceOrigin,
    pub completeness: AddgeneCompleteness,
    pub length: Option<usize>,
    /// Some catalogs provide only a file link; absence is never treated as sequence agreement.
    pub catalog_sequence_blake3: Option<String>,
}

#[derive(Clone)]
struct SequenceData {
    sequence: Option<Sequence>,
    download_url: Url,
}
/// Serializable summary excludes all presigned URLs and raw sequence data.
/// Keep this immutable value until selecting a sequence from the same response.
/// Catalog identity, category labels and response evidence cannot be rewritten.
/// Read-only getters preserve the exact parsed snapshot; owned choice copies are
/// presentation data and never replace the private import binding.
///
/// A caller cannot rewrite the plasmid identity:
/// ```compile_fail,E0616
/// fn forge(catalog: &mut collagene_sources::AddgeneCatalog) {
///     catalog.plasmid_id = 999;
/// }
/// ```
/// Nor can a caller rewrite its name:
/// ```compile_fail,E0616
/// fn forge(catalog: &mut collagene_sources::AddgeneCatalog) {
///     catalog.name = "Different plasmid".into();
/// }
/// ```
/// The response evidence hash is immutable:
/// ```compile_fail,E0616
/// fn forge(catalog: &mut collagene_sources::AddgeneCatalog) {
///     catalog.response_blake3 = "0".repeat(64);
/// }
/// ```
/// The ordered choice collection cannot be replaced:
/// ```compile_fail,E0616
/// fn forge(catalog: &mut collagene_sources::AddgeneCatalog) {
///     catalog.choices = Vec::new();
/// }
/// ```
/// A borrowed choice cannot be relabeled through its getter:
/// ```compile_fail,E0594
/// fn forge(catalog: &mut collagene_sources::AddgeneCatalog) {
///     catalog.choices()[0].origin = collagene_sources::AddgeneSequenceOrigin::Depositor;
/// }
/// ```
#[derive(Clone, Serialize)]
pub struct AddgeneCatalog {
    plasmid_id: u64,
    name: String,
    choices: Vec<AddgeneSequenceChoice>,
    /// Hash of the actual fetched catalog bytes, including its transient link representation.
    response_blake3: String,
    #[serde(skip)]
    records: std::sync::Arc<BTreeMap<u64, SequenceData>>,
}
impl fmt::Debug for AddgeneCatalog {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("AddgeneCatalog")
            .field("plasmid_id", &self.plasmid_id)
            .field("name", &self.name)
            .field("choices", &self.choices)
            .field("response_blake3", &self.response_blake3)
            .finish()
    }
}

pub fn addgene_catalog_url(plasmid_id: u64) -> Result<String, SourceError> {
    if plasmid_id == 0 {
        return Err(SourceError::InvalidSource);
    }
    Ok(format!(
        "https://{API_HOST}/catalog/plasmid-with-sequences/{plasmid_id}/"
    ))
}
fn string(value: &serde_json::Value, key: &str, max: usize) -> Result<String, SourceError> {
    let value = value
        .get(key)
        .and_then(|v| v.as_str())
        .ok_or(SourceError::AddgeneSchema)?;
    if value.len() > max
        || value
            .chars()
            .any(|c| c.is_control() && !matches!(c, '\n' | '\r' | '\t'))
    {
        return Err(SourceError::AddgeneSchema);
    }
    Ok(value.into())
}
#[derive(Deserialize)]
struct SequenceRecord {
    sequence_id: u64,
    sequence_description: String,
    length: Option<usize>,
    sequence: String,
    #[serde(default)]
    genbank_api_url: Option<String>,
    #[serde(default)]
    genbank_url: Option<String>,
}
/// Parse the documented `PlasmidWithSequences` response. Unknown metadata is ignored; unknown
/// sequence groups fail explicitly rather than hiding an additional record category.
pub fn parse_addgene_catalog(plasmid_id: u64, bytes: &[u8]) -> Result<AddgeneCatalog, SourceError> {
    if bytes.len() > MAX_IMPORT_BYTES {
        return Err(SourceError::TooLarge);
    }
    if plasmid_id == 0 {
        return Err(SourceError::InvalidSource);
    }
    let json: serde_json::Value =
        serde_json::from_slice(bytes).map_err(|_| SourceError::AddgeneSchema)?;
    if json.get("id").and_then(|v| v.as_u64()) != Some(plasmid_id) {
        return Err(SourceError::IdentifierMismatch);
    }
    let name = string(&json, "name", 1024)?;
    let groups = json
        .get("sequences")
        .and_then(|v| v.as_object())
        .ok_or(SourceError::AddgeneSchema)?;
    let definitions = [
        (
            "public_user_full_sequences",
            AddgeneSequenceOrigin::Depositor,
            AddgeneCompleteness::Full,
        ),
        (
            "public_addgene_full_sequences",
            AddgeneSequenceOrigin::Addgene,
            AddgeneCompleteness::Full,
        ),
        (
            "public_user_partial_sequences",
            AddgeneSequenceOrigin::Depositor,
            AddgeneCompleteness::Partial,
        ),
        (
            "public_addgene_partial_sequences",
            AddgeneSequenceOrigin::Addgene,
            AddgeneCompleteness::Partial,
        ),
    ];
    if groups
        .keys()
        .any(|key| !definitions.iter().any(|(name, _, _)| key == name))
    {
        return Err(SourceError::AddgeneSchema);
    }
    let mut catalog = AddgeneCatalog {
        plasmid_id,
        name,
        choices: vec![],
        response_blake3: blake3::hash(bytes).to_hex().to_string(),
        records: std::sync::Arc::new(BTreeMap::new()),
    };
    let mut seen = BTreeSet::new();
    for (key, origin, completeness) in definitions {
        let entries = groups
            .get(key)
            .and_then(|v| v.as_array())
            .ok_or(SourceError::AddgeneSchema)?;
        if entries.len() > MAX_CHOICES.saturating_sub(catalog.choices.len()) {
            return Err(SourceError::TooLarge);
        }
        for entry in entries {
            let record: SequenceRecord =
                serde_json::from_value(entry.clone()).map_err(|_| SourceError::AddgeneSchema)?;
            if record.sequence_id == 0
                || !seen.insert(record.sequence_id)
                || record.sequence_description.len() > 8192
                || record
                    .sequence_description
                    .chars()
                    .any(|c| c.is_control() && !matches!(c, '\n' | '\r' | '\t'))
                || record
                    .length
                    .is_some_and(|n| n == 0 || n > collagene_model::MAX_SEQUENCE_LEN)
            {
                return Err(SourceError::AddgeneSchema);
            }
            let sequence = if record.sequence.trim().is_empty() {
                None
            } else {
                Some(
                    Sequence::new(Alphabet::Dna, &record.sequence)
                        .map_err(|_| SourceError::InvalidSequenceData)?,
                )
            };
            if let Some(sequence) = &sequence {
                if record.length.is_some_and(|n| n != sequence.len()) {
                    return Err(SourceError::Truncated);
                }
            }
            // Current preferred key, with documented deprecated alias accepted only as a fallback.
            let raw_url = record
                .genbank_api_url
                .filter(|s| !s.trim().is_empty())
                .or_else(|| record.genbank_url.filter(|s| !s.trim().is_empty()))
                .ok_or(SourceError::AddgeneSchema)?;
            if raw_url.len() > 16384 {
                return Err(SourceError::TooLarge);
            }
            let download_url = Url::parse(&raw_url).map_err(|_| SourceError::AddgeneSchema)?;
            validate_public_url(&download_url)?;
            if download_url.host_str() == Some(API_HOST)
                && (download_url.path() != format!("/download/genbank/{}/", record.sequence_id)
                    || download_url.fragment().is_some())
            {
                return Err(SourceError::AddgeneSchema);
            }
            catalog.choices.push(AddgeneSequenceChoice {
                sequence_id: record.sequence_id,
                description: record.sequence_description,
                origin,
                completeness,
                length: record.length,
                catalog_sequence_blake3: sequence
                    .as_ref()
                    .map(|s| blake3::hash(s.as_str().as_bytes()).to_hex().to_string()),
            });
            std::sync::Arc::make_mut(&mut catalog.records).insert(
                record.sequence_id,
                SequenceData {
                    sequence,
                    download_url,
                },
            );
        }
    }
    Ok(catalog)
}
impl AddgeneCatalog {
    /// Exact plasmid identity checked against the fetched catalog response.
    pub fn plasmid_id(&self) -> u64 {
        self.plasmid_id
    }
    /// Provider-supplied name from the same immutable catalog snapshot.
    pub fn name(&self) -> &str {
        &self.name
    }
    /// Ordered, read-only provider choices. Clone individual records for display
    /// if needed; modifying a clone cannot change download/import provenance.
    pub fn choices(&self) -> &[AddgeneSequenceChoice] {
        &self.choices
    }
    /// BLAKE3 of the fetched response bytes, not proof of provider authenticity.
    pub fn response_blake3(&self) -> &str {
        &self.response_blake3
    }
    /// Transient presigned/provider download URL; do not save or log it.
    pub fn download_url(&self, sequence_id: u64) -> Result<String, SourceError> {
        self.records
            .get(&sequence_id)
            .map(|r| r.download_url.to_string())
            .ok_or(SourceError::AddgeneUnknownSequence)
    }
    /// Parse and check a selected record against the same immutable catalog response.
    /// API full/partial labels are provider evidence, not independent biological confirmation.
    pub fn import_download(
        &self,
        sequence_id: u64,
        bytes: &[u8],
        retrieved_at: &str,
    ) -> Result<ImportReport, SourceError> {
        let record = self
            .records
            .get(&sequence_id)
            .ok_or(SourceError::AddgeneUnknownSequence)?;
        let choice = self
            .choices
            .iter()
            .find(|c| c.sequence_id == sequence_id)
            .ok_or(SourceError::AddgeneUnknownSequence)?;
        if detect_format(bytes)? != DownloadFormat::GenBank {
            return Err(SourceError::UnexpectedContent);
        }
        let source = Source::Direct {
            url: format!("https://www.addgene.org/{}/", self.plasmid_id),
        };
        let mut report = import_download(&source, bytes, retrieved_at)?;
        if report.molecules.len() != 1 {
            return Err(SourceError::AmbiguousResponse {
                records: report.molecules.len(),
            });
        }
        let molecule = &mut report.molecules[0];
        if molecule.sequence.alphabet() != Alphabet::Dna
            || choice.length.is_some_and(|n| n != molecule.sequence.len())
            || record
                .sequence
                .as_ref()
                .is_some_and(|expected| *expected != molecule.sequence)
        {
            return Err(SourceError::AddgeneSequenceMismatch);
        }
        if choice.completeness == AddgeneCompleteness::Partial
            && molecule.topology == Topology::Circular
        {
            return Err(SourceError::AddgeneSequenceMismatch);
        }
        molecule
            .metadata
            .insert("source:provider".into(), "Addgene".into());
        molecule
            .metadata
            .insert("source:identifier".into(), self.plasmid_id.to_string());
        molecule
            .metadata
            .insert("source:sequence_id".into(), sequence_id.to_string());
        molecule.metadata.insert(
            "source:addgene_origin".into(),
            match choice.origin {
                AddgeneSequenceOrigin::Depositor => "depositor",
                AddgeneSequenceOrigin::Addgene => "addgene",
            }
            .into(),
        );
        molecule.metadata.insert(
            "source:addgene_completeness".into(),
            match choice.completeness {
                AddgeneCompleteness::Full => "full",
                AddgeneCompleteness::Partial => "partial",
            }
            .into(),
        );
        molecule.metadata.insert(
            "source:addgene_sequence_description".into(),
            choice.description.clone(),
        );
        molecule.metadata.insert(
            "source:addgene_catalog_blake3".into(),
            self.response_blake3.clone(),
        );
        // Keep only the supported, public facts selected by this adapter; never preserve URL tokens.
        molecule.metadata.insert(
            "source:addgene_catalog_evidence".into(),
            serde_json::to_string(
                &serde_json::json!({"plasmid_id":self.plasmid_id,"name":self.name,"choice":choice}),
            )
            .map_err(|_| SourceError::AddgeneSchema)?,
        );
        if record.sequence.is_none() {
            report.warnings.push("Addgene's catalog supplied no sequence bases for this record; the GenBank file could not be cross-checked against catalog bases.".into());
        }
        if choice.completeness == AddgeneCompleteness::Partial {
            report.warnings.push("This Addgene record is explicitly partial; it is not a complete plasmid reference.".into());
        }
        report.warnings.push("Addgene deposited/verified and full/partial categories are retained as provider assertions, not independent experimental verification.".into());
        collagene_core::protocol::save_document(molecule).map_err(|_| SourceError::TooLarge)?;
        Ok(report)
    }
}
