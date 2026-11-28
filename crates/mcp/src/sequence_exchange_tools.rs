//! Typed compatibility profiles, without host paths or restored chemistry authority.
use crate::DocumentRef;
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_core::{collagene_xml_compat, gcg_sequence};
use schemars::JsonSchema;
use serde::Deserialize;
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub enum ImportProfile {
    Gcg(gcg_sequence::ImportOptions),
    Fasta {
        /// None retains per-record declarations and defaults unmarked records to DNA.
        /// Some explicitly sets all records, refusing contradictory declarations.
        alphabet: Option<collagene_core::Alphabet>,
    },
    CollageneXml {
        policy: collagene_xml_compat::ReferencePolicy,
    },
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ImportQuery {
    pub text: String,
    pub profile: ImportProfile,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ReferenceExportQuery {
    pub document: DocumentRef,
    /// Default false returns a review instead of dropping unrepresentable fields.
    #[serde(default)]
    pub accept_omissions: bool,
}

pub const MAX_FILE_BYTES: usize = 8 * 1024 * 1024;
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub enum SequenceFileFormat {
    SnapGene,
    XdnaV0 {
        /// Fallback document name; the file's supported naming fields are retained.
        name: String,
    },
    Abif,
    Scf,
    NativeDocument,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct FileImportQuery {
    /// Exact file bytes, at most 8 MiB decoded. No host path or URL.
    pub data_base64: String,
    pub format: SequenceFileFormat,
}
pub fn import_file(
    args: FileImportQuery,
) -> Result<collagene_core::formats::ImportReport, collagene_core::Error> {
    use collagene_core::{Error, formats, legacy_formats, protocol, snapgene, traces};
    if args.data_base64.len() > MAX_FILE_BYTES.div_ceil(3) * 4 {
        return Err(Error::TooLarge);
    }
    let bytes = STANDARD
        .decode(args.data_base64)
        .map_err(|_| Error::Format("Sequence file data is not standard base64".into()))?;
    if bytes.len() > MAX_FILE_BYTES {
        return Err(Error::TooLarge);
    }
    match args.format {
        SequenceFileFormat::SnapGene => snapgene::import_snapgene(&bytes),
        SequenceFileFormat::XdnaV0 { name } => {
            if name.trim().is_empty() || name.len() > 256 || name.chars().any(char::is_control) {
                return Err(Error::Format("XDNA fallback name must be nonblank, at most 256 bytes and contain no controls".into()));
            }
            legacy_formats::import_xdna(&bytes, &name)
        }
        SequenceFileFormat::Abif | SequenceFileFormat::Scf => {
            let signature = match args.format {
                SequenceFileFormat::Abif => b"ABIF",
                _ => b".scf",
            };
            if !bytes.starts_with(signature) {
                return Err(Error::Format(
                    "Sequence file signature disagrees with the selected trace format".into(),
                ));
            }
            let trace = traces::import_trace(&bytes)?;
            Ok(formats::ImportReport {
                molecules: vec![trace.attached_molecule()?],
                warnings: trace.warnings,
            })
        }
        SequenceFileFormat::NativeDocument => Ok(formats::ImportReport {
            molecules: vec![protocol::load_document(&bytes)?],
            warnings: vec![],
        }),
    }
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct ScfExportQuery {
    pub document: DocumentRef,
    pub precision: collagene_core::scf_export::SamplePrecision,
    /// Documentary losses require review; this never permits invented signals/qualities.
    #[serde(default)]
    pub accept_omissions: bool,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct GcgExportQuery {
    pub document: DocumentRef,
    pub options: gcg_sequence::ExportOptions,
    /// Must explicitly acknowledge annotations/topology/history omitted by GCG.
    #[serde(default)]
    pub accept_omissions: bool,
}
