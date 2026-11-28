//! URL/source adapters, separate from the network-free scientific core.
//! Public resolution and provenance attachment also compile for browser WASM.
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_core::formats::sequence::{self as formats, ImportReport, MAX_IMPORT_BYTES};
use collagene_model::Molecule;
use std::{
    fmt,
    net::{IpAddr, Ipv4Addr, Ipv6Addr},
};
use url::{Host, Url};

mod addgene;
mod insdc;
pub use addgene::{
    AddgeneCatalog, AddgeneCompleteness, AddgeneSequenceChoice, AddgeneSequenceOrigin,
    addgene_catalog_url, parse_addgene_catalog,
};

#[cfg(not(target_arch = "wasm32"))]
mod native;
#[cfg(not(target_arch = "wasm32"))]
pub use native::SourceClient;

#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub enum SourceError {
    #[error(
        "Enter an NCBI accession/version, ena:/ddbj: accession, Addgene ID/page URL, or HTTPS sequence-file URL"
    )]
    InvalidSource,
    #[error(
        "Only public HTTPS URLs on the standard TLS port are supported; URL user credentials and local/private destinations are forbidden"
    )]
    DisallowedUrl,
    #[error(
        "Search/range/partial/reverse/collection selectors are not single complete-record imports; supply one accession.version"
    )]
    UnsupportedSelection,
    #[error(
        "Addgene plasmid {plasmid_id} requires approved API/MCP sequence access. Configure an approved Addgene API token, or provide an official GenBank download link or export instead"
    )]
    AddgeneAccessRequired { plasmid_id: u64 },
    #[error(
        "Addgene plasmid {plasmid_id} has {sequences} sequence choices; select a sequence_id explicitly"
    )]
    AddgeneSelectionRequired { plasmid_id: u64, sequences: usize },
    #[error("Addgene returned no available sequence records for plasmid {plasmid_id}")]
    AddgeneNoSequences { plasmid_id: u64 },
    #[error("The selected Addgene sequence ID is absent from this plasmid catalog response")]
    AddgeneUnknownSequence,
    #[error("Addgene catalog response does not match the supported documented schema")]
    AddgeneSchema,
    #[error("Addgene catalog sequence and downloaded GenBank sequence disagree")]
    AddgeneSequenceMismatch,
    #[error("Sequence download exceeded the 32 MiB limit")]
    TooLarge,
    #[error("Sequence download timed out")]
    Timeout,
    #[error(
        "Sequence download failed; no URL credentials or response body are included in this error"
    )]
    DownloadFailed,
    #[error("Sequence source returned HTTP {status}; retry after seconds: {retry_after_seconds:?}")]
    HttpStatus {
        status: u16,
        retry_after_seconds: Option<u64>,
    },
    #[error("Sequence source rate limit reached; retry after seconds: {retry_after_seconds:?}")]
    RateLimited { retry_after_seconds: Option<u64> },
    #[error("Sequence source returned a login/HTML/error page or an unsupported file type")]
    UnexpectedContent,
    #[error("Sequence download was incomplete or disagreed with its declared length")]
    Truncated,
    #[error("Downloaded sequence data did not pass the shared format parser")]
    InvalidSequenceData,
    #[error("Sequence provider did not return a valid accession.version")]
    MissingAccessionVersion,
    #[error("Returned accession.version differs from the requested identifier")]
    IdentifierMismatch,
    #[error(
        "A single provider record request returned {records} records; an explicit selection is required"
    )]
    AmbiguousResponse { records: usize },
    #[error("Too many redirects or a redirect attempted to leave the permitted source boundary")]
    RedirectRejected,
    #[error("Invalid source-client configuration")]
    InvalidConfig,
}

/// Transient source locator. Debug output deliberately omits signed query strings.
#[derive(Clone, PartialEq, Eq)]
pub enum Source {
    Ncbi { accession: String },
    Ena { accession: String },
    Ddbj { accession: String },
    Addgene { plasmid_id: u64 },
    Direct { url: String },
}
impl fmt::Debug for Source {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Ena { accession: value } | Self::Ddbj { accession: value } => f
                .debug_struct(if matches!(self, Self::Ena { .. }) {
                    "Ena"
                } else {
                    "Ddbj"
                })
                .field(
                    "accession",
                    &insdc::sequence_accession(value).unwrap_or_else(|| "<invalid>".into()),
                )
                .finish(),
            Self::Ncbi { accession: value } => f
                .debug_struct("Ncbi")
                .field(
                    "accession",
                    &accession(value, true).unwrap_or_else(|| "<invalid>".into()),
                )
                .finish(),
            Self::Addgene { plasmid_id } => f
                .debug_struct("Addgene")
                .field("plasmid_id", plasmid_id)
                .finish(),
            Self::Direct { .. } => f
                .debug_struct("Direct")
                .field("canonical_url", &self.canonical_url())
                .finish(),
        }
    }
}

/// Credentials are private, neither Debug nor serializable, and scoped to provider origins.
#[derive(Clone)]
pub struct Config {
    ncbi_api_key: Option<String>,
    contact_email: Option<String>,
    addgene_api_token: Option<String>,
}
impl Default for Config {
    fn default() -> Self {
        Self {
            ncbi_api_key: None,
            contact_email: None,
            addgene_api_token: None,
        }
    }
}
impl Config {
    pub fn with_addgene_api_token(mut self, token: impl Into<String>) -> Result<Self, SourceError> {
        let token = token.into();
        if token.is_empty() || token.len() > 4096 || !token.bytes().all(|b| b.is_ascii_graphic()) {
            return Err(SourceError::InvalidConfig);
        }
        self.addgene_api_token = Some(token);
        Ok(self)
    }
    pub fn with_ncbi_api_key(mut self, key: impl Into<String>) -> Result<Self, SourceError> {
        let key = key.into();
        if key.is_empty()
            || key.len() > 256
            || !key
                .bytes()
                .all(|b| b.is_ascii_alphanumeric() || b == b'-' || b == b'_')
        {
            return Err(SourceError::InvalidConfig);
        }
        self.ncbi_api_key = Some(key);
        Ok(self)
    }
    pub fn with_contact_email(mut self, email: impl Into<String>) -> Result<Self, SourceError> {
        let email = email.into();
        if email.len() > 254
            || !email.contains('@')
            || email
                .bytes()
                .any(|b| b.is_ascii_whitespace() || b.is_ascii_control())
        {
            return Err(SourceError::InvalidConfig);
        }
        self.contact_email = Some(email);
        Ok(self)
    }
    #[cfg(not(target_arch = "wasm32"))]
    pub fn from_env() -> Result<Self, SourceError> {
        let mut config = Self::default();
        if let Ok(key) = std::env::var("COLLAGENE_NCBI_API_KEY") {
            config = config.with_ncbi_api_key(key)?;
        }
        if let Ok(email) = std::env::var("COLLAGENE_NCBI_EMAIL") {
            config = config.with_contact_email(email)?;
        }
        if let Ok(token) = std::env::var("COLLAGENE_ADDGENE_API_TOKEN") {
            config = config.with_addgene_api_token(token)?;
        }
        Ok(config)
    }
}
fn accession(value: &str, allow_numeric: bool) -> Option<String> {
    if value.is_empty() || value.len() > 64 || !value.is_ascii() {
        return None;
    }
    if value.bytes().all(|b| b.is_ascii_digit()) {
        return allow_numeric
            .then(|| value.parse::<u64>().ok())
            .flatten()
            .filter(|n| *n > 0)
            .map(|_| value.into());
    }
    let mut parts = value.split('.');
    let base = parts.next()?;
    if let Some(version) = parts.next() {
        if version.is_empty()
            || !version.bytes().all(|b| b.is_ascii_digit())
            || version.parse::<u32>().ok()? == 0
            || parts.next().is_some()
        {
            return None;
        }
    }
    if !base.as_bytes().first()?.is_ascii_alphabetic()
        || !base.bytes().all(|b| b.is_ascii_alphanumeric() || b == b'_')
        || !base.bytes().any(|b| b.is_ascii_digit())
    {
        return None;
    }
    Some(value.to_ascii_uppercase())
}
fn plasmid_id(value: &str) -> Option<u64> {
    if value.bytes().all(|b| b.is_ascii_digit()) {
        value.parse().ok().filter(|n| *n > 0)
    } else {
        None
    }
}
fn ipv4_public(ip: Ipv4Addr) -> bool {
    let [a, b, c, _] = ip.octets();
    !ip.is_private()
        && !ip.is_loopback()
        && !ip.is_link_local()
        && !ip.is_documentation()
        && !ip.is_unspecified()
        && !ip.is_multicast()
        && !ip.is_broadcast()
        && a != 0
        && a < 224
        && !(a == 100 && (64..=127).contains(&b))
        && !(a == 192 && b == 0 && c == 0)
        && !(a == 198 && (b == 18 || b == 19))
        && !(a == 192 && b == 88 && c == 99)
}
fn ipv6_public(ip: Ipv6Addr) -> bool {
    let s = ip.segments();
    // Conservative global-unicast-only policy, excluding transition/documentation/special ranges.
    s[0] & 0xe000 == 0x2000
        && s[0] != 0x2002
        && s[0] != 0x3fff
        && !(s[0] == 0x2001 && (s[1] < 0x0200 || s[1] == 0x0db8))
}
pub fn is_public_ip(ip: IpAddr) -> bool {
    match ip {
        IpAddr::V4(v) => ipv4_public(v),
        IpAddr::V6(v) => ipv6_public(v),
    }
}
pub fn validate_public_url(url: &Url) -> Result<(), SourceError> {
    if url.scheme() != "https"
        || !url.username().is_empty()
        || url.password().is_some()
        || url.port().is_some_and(|p| p != 443)
    {
        return Err(SourceError::DisallowedUrl);
    }
    match url.host().ok_or(SourceError::DisallowedUrl)? {
        Host::Ipv4(ip) if !ipv4_public(ip) => return Err(SourceError::DisallowedUrl),
        Host::Ipv6(ip) if !ipv6_public(ip) => return Err(SourceError::DisallowedUrl),
        Host::Domain(host) => {
            let host = host.trim_end_matches('.').to_ascii_lowercase();
            if !host.contains('.')
                || host == "localhost"
                || [
                    ".localhost",
                    ".local",
                    ".internal",
                    ".home",
                    ".test",
                    ".invalid",
                ]
                .iter()
                .any(|suffix| host.ends_with(suffix))
            {
                return Err(SourceError::DisallowedUrl);
            }
        }
        _ => {}
    }
    Ok(())
}
fn sanitized(url: &Url) -> String {
    let mut clean = url.clone();
    clean.set_query(None);
    clean.set_fragment(None);
    let _ = clean.set_username("");
    let _ = clean.set_password(None);
    clean.to_string()
}

pub fn parse_source(input: &str) -> Result<Source, SourceError> {
    let input = input.trim();
    if input.is_empty() || input.len() > 8192 {
        return Err(SourceError::InvalidSource);
    }
    for (prefix, ena) in [
        ("ena:", true),
        ("ENA:", true),
        ("ddbj:", false),
        ("DDBJ:", false),
    ] {
        if let Some(id) = input.strip_prefix(prefix) {
            let accession =
                insdc::sequence_accession(id.trim()).ok_or(SourceError::InvalidSource)?;
            return Ok(if ena {
                Source::Ena { accession }
            } else {
                Source::Ddbj { accession }
            });
        }
    }
    if let Some(id) = input
        .strip_prefix("addgene:")
        .or_else(|| input.strip_prefix("Addgene:"))
    {
        return plasmid_id(id.trim())
            .map(|plasmid_id| Source::Addgene { plasmid_id })
            .ok_or(SourceError::InvalidSource);
    }
    if let Some(id) = input
        .strip_prefix("ncbi:")
        .or_else(|| input.strip_prefix("NCBI:"))
    {
        return accession(id.trim(), true)
            .map(|accession| Source::Ncbi { accession })
            .ok_or(SourceError::InvalidSource);
    }
    if let Some(plasmid_id) = plasmid_id(input.trim_start_matches('#')) {
        return Ok(Source::Addgene { plasmid_id });
    }
    if !input.contains("://") {
        return accession(input, false)
            .map(|accession| Source::Ncbi { accession })
            .ok_or(SourceError::InvalidSource);
    }
    let mut url = Url::parse(input).map_err(|_| SourceError::InvalidSource)?;
    validate_public_url(&url)?;
    url.set_fragment(None);
    if let Some(result) = insdc::resolve(&url) {
        return result;
    }
    let host = url
        .host_str()
        .ok_or(SourceError::InvalidSource)?
        .to_ascii_lowercase();
    let path: Vec<_> = url
        .path_segments()
        .ok_or(SourceError::InvalidSource)?
        .filter(|s| !s.is_empty())
        .collect();
    if matches!(host.as_str(), "www.addgene.org" | "addgene.org")
        && (path.len() == 1 || path.len() == 2 && path[1] == "sequences")
    {
        if let Some(plasmid_id) = path.first().and_then(|s| plasmid_id(s)) {
            return Ok(Source::Addgene { plasmid_id });
        }
    }
    if matches!(
        host.as_str(),
        "www.ncbi.nlm.nih.gov" | "ncbi.nlm.nih.gov" | "eutils.ncbi.nlm.nih.gov"
    ) {
        let queries: Vec<_> = url.query_pairs().collect();
        if queries.iter().any(|(k, _)| {
            matches!(
                k.as_ref(),
                "seq_start"
                    | "seq_stop"
                    | "from"
                    | "to"
                    | "strand"
                    | "term"
                    | "query_key"
                    | "WebEnv"
            )
        }) {
            return Err(SourceError::UnsupportedSelection);
        }
        let id = if path.len() == 2 && matches!(path[0], "nuccore" | "nucleotide") {
            Some(path[1].to_string())
        } else if url.path() == "/entrez/eutils/efetch.fcgi" || url.path() == "/sviewer/viewer.fcgi"
        {
            if queries
                .iter()
                .any(|(k, v)| k == "db" && !matches!(v.as_ref(), "nuccore" | "nucleotide"))
            {
                return Err(SourceError::UnsupportedSelection);
            }
            let ids: Vec<_> = queries
                .iter()
                .filter(|(k, _)| k == "id")
                .map(|(_, v)| v.to_string())
                .collect();
            if ids.len() != 1 {
                return Err(SourceError::UnsupportedSelection);
            }
            ids.into_iter().next()
        } else {
            None
        };
        let id = id.ok_or(SourceError::UnsupportedSelection)?;
        return accession(&id, true)
            .map(|accession| Source::Ncbi { accession })
            .ok_or(SourceError::UnsupportedSelection);
    }
    Ok(Source::Direct {
        url: url.to_string(),
    })
}
impl Source {
    pub fn canonical_url(&self) -> String {
        match self {
            Self::Ncbi { accession } => format!("https://www.ncbi.nlm.nih.gov/nuccore/{accession}"),
            Self::Ena { accession } => {
                format!("https://www.ebi.ac.uk/ena/browser/view/{accession}")
            }
            Self::Ddbj { accession } => {
                format!("https://getentry.ddbj.nig.ac.jp/getentry/na/{accession}")
            }
            Self::Addgene { plasmid_id } => format!("https://www.addgene.org/{plasmid_id}/"),
            Self::Direct { url } => Url::parse(url).map(|u| sanitized(&u)).unwrap_or_default(),
        }
    }
    /// Transient fetch target, possibly containing a configured API key. Never persist/log this URL.
    pub fn download_url(&self) -> Result<String, SourceError> {
        Ok(self
            .download_url_with_config(&Config::default())?
            .to_string())
    }
    fn download_url_with_config(&self, config: &Config) -> Result<Url, SourceError> {
        if let Self::Ncbi { accession: value } = self {
            if accession(value, true).as_deref() != Some(value.as_str()) {
                return Err(SourceError::InvalidSource);
            }
        }
        match self {
            Self::Ena { .. } | Self::Ddbj { .. } => insdc::download(self),
            Self::Ncbi { accession } => {
                let mut url =
                    Url::parse("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi")
                        .unwrap();
                {
                    let mut q = url.query_pairs_mut();
                    q.append_pair("db", "nuccore")
                        .append_pair("id", accession)
                        .append_pair("rettype", "gbwithparts")
                        .append_pair("retmode", "text")
                        .append_pair("complexity", "1")
                        .append_pair("tool", "collagene2");
                    if let Some(email) = &config.contact_email {
                        q.append_pair("email", email);
                    }
                    if let Some(key) = &config.ncbi_api_key {
                        q.append_pair("api_key", key);
                    }
                }
                Ok(url)
            }
            Self::Addgene { plasmid_id } => Err(SourceError::AddgeneAccessRequired {
                plasmid_id: *plasmid_id,
            }),
            Self::Direct { url } => {
                let u = Url::parse(url).map_err(|_| SourceError::InvalidSource)?;
                validate_public_url(&u)?;
                Ok(u)
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DownloadFormat {
    GenBank,
    Fasta,
    Embl,
    SnapGene,
    Abif,
    Scf,
    Collagene,
}
impl DownloadFormat {
    fn name(self) -> &'static str {
        match self {
            Self::GenBank => "GenBank",
            Self::Fasta => "FASTA",
            Self::Embl => "EMBL",
            Self::SnapGene => "SnapGene DNA",
            Self::Abif => "ABIF",
            Self::Scf => "SCF",
            Self::Collagene => "Collagene 2",
        }
    }
}
pub fn detect_format(bytes: &[u8]) -> Result<DownloadFormat, SourceError> {
    if bytes.len() > MAX_IMPORT_BYTES {
        return Err(SourceError::TooLarge);
    }
    if bytes.starts_with(b"ABIF") {
        return Ok(DownloadFormat::Abif);
    }
    if bytes.starts_with(b".scf") {
        return Ok(DownloadFormat::Scf);
    }
    if bytes.starts_with(b"CG2\0") {
        return Ok(DownloadFormat::Collagene);
    }
    if bytes.len() >= 13 && &bytes[5..13] == b"SnapGene" {
        return Ok(DownloadFormat::SnapGene);
    }
    let text = std::str::from_utf8(bytes)
        .map_err(|_| SourceError::UnexpectedContent)?
        .trim_start_matches('\u{feff}')
        .trim_start();
    if text.starts_with("LOCUS ") {
        Ok(DownloadFormat::GenBank)
    } else if text.starts_with("ID ") {
        Ok(DownloadFormat::Embl)
    } else if text.starts_with('>') {
        Ok(DownloadFormat::Fasta)
    } else {
        Err(SourceError::UnexpectedContent)
    }
}
fn returned_accession(m: &Molecule) -> Option<String> {
    m.original_genbank
        .as_ref()
        .and_then(|raw| raw.lines().find(|line| line.starts_with("VERSION ")))
        .and_then(|line| line.split_ascii_whitespace().nth(1))
        .and_then(|id| accession(id, false))
        .filter(|id| id.contains('.'))
}
fn validate_genbank_lengths(molecules: &[Molecule]) -> Result<(), SourceError> {
    for m in molecules {
        let raw = m
            .original_genbank
            .as_ref()
            .ok_or(SourceError::InvalidSequenceData)?;
        if !raw.trim_end().ends_with("//") {
            return Err(SourceError::Truncated);
        }
        let declared = raw
            .lines()
            .find(|l| l.starts_with("LOCUS "))
            .and_then(|l| l.split_ascii_whitespace().nth(2))
            .and_then(|n| n.parse::<usize>().ok())
            .ok_or(SourceError::InvalidSequenceData)?;
        if declared != m.sequence.len() {
            return Err(SourceError::Truncated);
        }
    }
    Ok(())
}
fn fasta_records(bytes: &[u8]) -> Vec<&[u8]> {
    let mut starts = vec![];
    let mut offset = 0;
    for line in bytes.split_inclusive(|b| *b == b'\n') {
        let line = if offset == 0 {
            line.strip_prefix(&[0xef, 0xbb, 0xbf]).unwrap_or(line)
        } else {
            line
        };
        if line.trim_ascii_start().starts_with(b">") {
            starts.push(offset);
        }
        offset += if offset == 0 && bytes.starts_with(&[0xef, 0xbb, 0xbf]) {
            line.len() + 3
        } else {
            line.len()
        };
    }
    if let Some(first) = starts.first_mut() {
        *first = 0;
    }
    starts
        .iter()
        .enumerate()
        .map(|(i, start)| &bytes[*start..starts.get(i + 1).copied().unwrap_or(bytes.len())])
        .collect()
}

/// Apply the shared parser and attach non-secret provenance after a completed download.
/// Browser callers perform their own fetch/CORS handling and pass an actual UTC retrieval time.
pub fn annotate_import(
    report: &mut ImportReport,
    source: &Source,
    bytes: &[u8],
    retrieved_at: &str,
) -> Result<(), SourceError> {
    if retrieved_at.is_empty() || retrieved_at.len() > 64 || !retrieved_at.is_ascii() {
        return Err(SourceError::InvalidConfig);
    }
    source.download_url()?;
    let retrieved_at = chrono::DateTime::parse_from_rfc3339(retrieved_at)
        .map_err(|_| SourceError::InvalidConfig)?
        .with_timezone(&chrono::Utc)
        .to_rfc3339_opts(chrono::SecondsFormat::AutoSi, true);
    let format = detect_format(bytes)?;
    if matches!(source, Source::Addgene { .. }) {
        return Err(SourceError::UnexpectedContent);
    }
    if matches!(source, Source::Ncbi { .. } | Source::Ddbj { .. })
        && format != DownloadFormat::GenBank
        || matches!(source, Source::Ena { .. }) && format != DownloadFormat::Embl
    {
        return Err(SourceError::UnexpectedContent);
    }
    if format == DownloadFormat::GenBank {
        validate_genbank_lengths(&report.molecules)?;
    }
    if matches!(
        source,
        Source::Ncbi { .. } | Source::Ena { .. } | Source::Ddbj { .. }
    ) && report.molecules.len() != 1
    {
        return Err(SourceError::AmbiguousResponse {
            records: report.molecules.len(),
        });
    }
    // The portable annotation entry point cannot attach provider authority to an
    // edited or unrelated caller-supplied document using only an old raw record.
    if matches!(source, Source::Ena { .. } | Source::Ddbj { .. }) {
        let parsed = formats::import_bytes(bytes).map_err(|_| SourceError::InvalidSequenceData)?;
        if parsed.molecules != report.molecules {
            return Err(SourceError::InvalidSequenceData);
        }
    }
    let content_hash = blake3::hash(bytes).to_hex().to_string();
    let fasta = if format == DownloadFormat::Fasta {
        fasta_records(bytes)
    } else {
        vec![]
    };
    for (index, m) in report.molecules.iter_mut().enumerate() {
        let actual = match format {
            DownloadFormat::GenBank => returned_accession(m),
            DownloadFormat::Embl if matches!(source, Source::Ena { .. }) => {
                Some(insdc::actual_embl(m)?)
            }
            DownloadFormat::Embl => m
                .metadata
                .get("embl.sequence_version")
                .filter(|v| !v.is_empty())
                .and_then(|v| accession(&format!("{}.{v}", m.name), false)),
            _ => None,
        };
        m.metadata.remove("source:accession_version");
        let canonical = match source {
            Source::Ncbi {
                accession: requested,
            }
            | Source::Ena {
                accession: requested,
            }
            | Source::Ddbj {
                accession: requested,
            } => {
                let actual = actual
                    .as_ref()
                    .ok_or(SourceError::MissingAccessionVersion)?;
                if requested.bytes().any(|b| b.is_ascii_alphabetic())
                    && (requested.contains('.') && requested != actual
                        || !requested.contains('.')
                            && requested != actual.split('.').next().unwrap_or_default())
                {
                    return Err(SourceError::IdentifierMismatch);
                }
                if matches!(source, Source::Ddbj { .. }) {
                    insdc::validate_genbank_identity(m, actual)?;
                }
                match source {
                    Source::Ena { .. } => Source::Ena {
                        accession: actual.clone(),
                    }
                    .canonical_url(),
                    Source::Ddbj { .. } => Source::Ddbj {
                        accession: actual.clone(),
                    }
                    .canonical_url(),
                    _ => format!("https://www.ncbi.nlm.nih.gov/nuccore/{actual}"),
                }
            }
            _ => source.canonical_url(),
        };
        m.metadata.insert("source:canonical_url".into(), canonical);
        m.metadata.insert(
            "source:provider".into(),
            match source {
                Source::Ncbi { .. } => "NCBI nuccore",
                Source::Ena { .. } => "ENA",
                Source::Ddbj { .. } => "DDBJ",
                Source::Addgene { .. } => "Addgene",
                Source::Direct { .. } => "Direct HTTPS",
            }
            .into(),
        );
        m.metadata
            .insert("source:retrieved_at".into(), retrieved_at.clone());
        m.metadata
            .insert("source:content_blake3".into(), content_hash.clone());
        m.metadata
            .insert("source:format".into(), format.name().into());
        m.metadata
            .insert("source:record_index".into(), index.to_string());
        if let Source::Ncbi { accession } | Source::Ena { accession } | Source::Ddbj { accession } =
            source
        {
            m.metadata
                .insert("source:identifier".into(), accession.clone());
        } else {
            m.metadata
                .insert("source:identifier".into(), content_hash.clone());
        }
        if let Some(actual) = actual {
            m.metadata.insert("source:accession_version".into(), actual);
        }
        if format == DownloadFormat::Fasta {
            let record = fasta
                .get(index)
                .copied()
                .ok_or(SourceError::InvalidSequenceData)?;
            m.metadata.insert(
                "source:original_fasta_base64".into(),
                STANDARD.encode(record),
            );
            m.metadata.insert(
                "source:original_fasta_blake3".into(),
                blake3::hash(record).to_hex().to_string(),
            );
        }
        // Imported evidence and provenance must still fit a native document payload.
        collagene_core::protocol::save_document(m).map_err(|_| SourceError::TooLarge)?;
    }
    if format == DownloadFormat::Fasta {
        report.warnings.push("FASTA does not declare sequence lengths; transport completion does not establish biological completeness. Original record bytes were retained.".into());
    }
    report.warnings.push("URL-import provenance identifies fetched data; it does not independently verify a complete plasmid, annotations or experimental function. Re-import creates documents rather than overwriting an existing design.".into());
    Ok(())
}
pub fn import_download(
    source: &Source,
    bytes: &[u8],
    retrieved_at: &str,
) -> Result<ImportReport, SourceError> {
    detect_format(bytes)?;
    let mut report = formats::import_bytes(bytes).map_err(|_| SourceError::InvalidSequenceData)?;
    annotate_import(&mut report, source, bytes, retrieved_at)?;
    Ok(report)
}
