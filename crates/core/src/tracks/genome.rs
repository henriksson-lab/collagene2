//! Immutable large-reference blocks and streaming FASTA ingestion, independent of Molecule.
use crate::sequence::analysis;
use collagene_model::{Error, Topology};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
pub const BLOCK_BASES: usize = 65_536;
pub const HEADER: [u8; 8] = *b"CGREF\0\x01\0";
pub const TAIL_BYTES: usize = 48;
pub const MAX_MANIFEST_BYTES: usize = 16 * 1024 * 1024;
const TAIL_MAGIC: [u8; 8] = *b"CGREND1\0";
const MAX_CONTIGS: usize = 10_000;
const MAX_BLOCKS: usize = 250_000;
macro_rules! dto{($($x:tt)*)=>{#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($x)*};}
dto! {pub struct Provenance{pub source:String,pub identifier:Option<String>,pub retrieved_at:Option<String>,pub content_hash:Option<String>}}
dto! {pub struct ExpectedComponent{pub id:String,pub source_id:String,pub accession:Option<String>,pub length:u64,pub sequence_hash:Option<String>,pub topology:Topology}}
dto! {pub struct BuildInfo{pub name:String,pub assembly_accession:Option<String>,pub version:Option<String>,pub provenance:Vec<Provenance>,pub expected_components:Vec<ExpectedComponent>}}
impl Default for BuildInfo {
    fn default() -> Self {
        Self {
            name: "Imported reference collection".into(),
            assembly_accession: None,
            version: None,
            provenance: vec![],
            expected_components: vec![],
        }
    }
}
dto! {pub struct IngestOptions{pub max_source_bytes:u64,pub max_total_bases:u64,pub max_contig_bases:u64,pub max_contigs:usize,pub max_blocks:usize,pub max_push_bytes:usize}}
impl Default for IngestOptions {
    fn default() -> Self {
        Self {
            max_source_bytes: 32 * 1024 * 1024 * 1024,
            max_total_bases: 16_000_000_000,
            max_contig_bases: 10_000_000_000,
            max_contigs: 10_000,
            max_blocks: 250_000,
            max_push_bytes: 1024 * 1024,
        }
    }
}
dto! {pub struct BlockDescriptor{pub index:u64,pub start:u64,pub count:u32,pub offset:u64,pub encoded_hash:String}}
dto! {pub struct ContigRecord{pub id:String,pub source_id:String,pub accession:Option<String>,pub description:String,pub length:u64,pub topology:Topology,pub sequence_hash:String,pub blocks:Vec<BlockDescriptor>}}
dto! {pub struct GenomeManifest{pub schema_version:u32,pub build:BuildInfo,pub contigs:Vec<ContigRecord>,pub content_revision:String,pub raw_source_hash:String,pub source_bytes:u64,pub total_bases:u64,pub missing_expected:Vec<String>,pub expected_complete:bool}}
dto! {pub struct EncodedBlock{pub contig_id:String,pub descriptor:BlockDescriptor,pub encoded:Vec<u8>}}
dto! {pub struct FinishedReference{pub manifest:GenomeManifest,pub manifest_hash:String,pub manifest_offset:u64,pub file_length:u64}}
dto! {pub struct TailInfo{pub manifest_length:u64,pub manifest_hash:String,pub manifest_offset:u64,pub file_length:u64}}
dto! {pub struct GenomicInterval{pub contig_id:String,pub start:u64,pub end:u64}}
dto! {pub struct GlobalPosition{pub content_revision:String,pub contig_id:String,pub position:u64}}
dto! {pub struct SemanticReceipt{pub manifest_hash:String,pub content_revision:String,pub contigs_verified:usize,pub total_bases:u64}}
pub const WORKSPACE_ATTACHMENT_KIND: &str = "genome-reference-link-v1";
dto! {pub struct ReferenceLink{pub schema_version:u32,pub role:String,pub manifest_hash:String,pub content_revision:String,pub name:String,pub assembly_accession:Option<String>,pub version:Option<String>,pub contig_count:usize,pub total_bases:u64,pub file_length:u64,pub cache_hint:Option<String>}}
fn bad(s: &str) -> Error {
    Error::Format(format!("genome reference: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_hexdigit() && !b.is_ascii_uppercase())
}
fn text(s: &str, max: usize) -> bool {
    !s.is_empty() && s.len() <= max && !s.chars().any(char::is_control)
}
fn id(s: &str) -> bool {
    text(s, 256) && !s.chars().any(char::is_whitespace)
}
fn optional(s: &Option<String>, max: usize) -> bool {
    s.as_ref().is_none_or(|s| text(s, max))
}
struct Payload(Vec<u8>);
impl postcard::ser_flavors::Flavor for Payload {
    type Output = Vec<u8>;
    fn try_push(&mut self, b: u8) -> postcard::Result<()> {
        if self.0.len() >= MAX_MANIFEST_BYTES {
            return Err(postcard::Error::SerializeBufferFull);
        }
        self.0.push(b);
        Ok(())
    }
    fn try_extend(&mut self, b: &[u8]) -> postcard::Result<()> {
        if b.len() > MAX_MANIFEST_BYTES.saturating_sub(self.0.len()) {
            return Err(postcard::Error::SerializeBufferFull);
        }
        self.0.extend_from_slice(b);
        Ok(())
    }
    fn finalize(self) -> postcard::Result<Vec<u8>> {
        Ok(self.0)
    }
}
fn encode(v: &impl Serialize) -> Result<Vec<u8>, Error> {
    postcard::serialize_with_flavor(v, Payload(Vec::new())).map_err(|_| Error::TooLarge)
}
impl IngestOptions {
    pub fn validate(&self) -> Result<(), Error> {
        let bases = MAX_BLOCKS as u64 * BLOCK_BASES as u64;
        if self.max_source_bytes == 0
            || self.max_source_bytes > 128 * 1024 * 1024 * 1024
            || self.max_total_bases == 0
            || self.max_total_bases > bases
            || self.max_contig_bases == 0
            || self.max_contig_bases > bases
            || self.max_contigs == 0
            || self.max_contigs > MAX_CONTIGS
            || self.max_blocks == 0
            || self.max_blocks > MAX_BLOCKS
            || self.max_push_bytes == 0
            || self.max_push_bytes > 1024 * 1024
        {
            return Err(bad("invalid bounded ingestion options"));
        }
        Ok(())
    }
}
impl BuildInfo {
    pub fn validate(&self) -> Result<(), Error> {
        if !text(&self.name, 4096)
            || !optional(&self.assembly_accession, 256)
            || !optional(&self.version, 256)
            || self.provenance.len() > 128
            || self.expected_components.len() > MAX_CONTIGS
        {
            return Err(bad("invalid build metadata bounds"));
        }
        for p in &self.provenance {
            if !text(&p.source, 4096)
                || !optional(&p.identifier, 512)
                || !optional(&p.retrieved_at, 128)
                || p.content_hash.as_ref().is_some_and(|h| !hex(h))
            {
                return Err(bad("invalid source provenance"));
            }
        }
        let mut ids = BTreeSet::new();
        let mut sources = BTreeSet::new();
        for c in &self.expected_components {
            if !id(&c.id)
                || !id(&c.source_id)
                || !ids.insert(&c.id)
                || !sources.insert(&c.source_id)
                || !optional(&c.accession, 256)
                || c.length == 0
                || c.length > MAX_BLOCKS as u64 * BLOCK_BASES as u64
                || c.sequence_hash.as_ref().is_some_and(|h| !hex(h))
            {
                return Err(bad(
                    "expected component identity/length/hash invalid or duplicated",
                ));
            }
        }
        encode(self)?;
        Ok(())
    }
}
fn content_revision(contigs: &[ContigRecord]) -> Result<String, Error> {
    let mut roles = contigs
        .iter()
        .map(|c| (&c.id, &c.accession, c.length, c.topology, &c.sequence_hash))
        .collect::<Vec<_>>();
    roles.sort_by(|a, b| a.0.cmp(b.0));
    Ok(
        blake3::hash(&encode(&("cgref-canonical-content-v1", roles))?)
            .to_hex()
            .to_string(),
    )
}
fn missing(build: &BuildInfo, contigs: &[ContigRecord]) -> Result<Vec<String>, Error> {
    let mut absent = vec![];
    for e in &build.expected_components {
        if let Some(c) = contigs.iter().find(|c| c.source_id == e.source_id) {
            if c.id != e.id
                || c.accession != e.accession
                || c.length != e.length
                || c.topology != e.topology
                || e.sequence_hash
                    .as_ref()
                    .is_some_and(|h| h != &c.sequence_hash)
            {
                return Err(bad(
                    "observed component contradicts expected identity, length, topology or canonical hash",
                ));
            }
        } else {
            absent.push(e.id.clone());
        }
    }
    Ok(absent)
}
impl GenomeManifest {
    /// Bounded layout/metadata validation. Does NOT attest canonical DNA content.
    pub fn validate(&self) -> Result<(), Error> {
        self.build.validate()?;
        if self.schema_version != 1
            || self.contigs.is_empty()
            || self.contigs.len() > MAX_CONTIGS
            || !hex(&self.content_revision)
            || !hex(&self.raw_source_hash)
            || self.source_bytes == 0
            || self.source_bytes > 128 * 1024 * 1024 * 1024
            || self.missing_expected.len() > MAX_CONTIGS
        {
            return Err(bad("invalid manifest header"));
        }
        let mut roles = BTreeSet::new();
        let mut sources = BTreeSet::new();
        let mut total = 0u64;
        let mut offset = HEADER.len() as u64;
        let mut blocks = 0usize;
        for c in &self.contigs {
            if !id(&c.id)
                || !id(&c.source_id)
                || !roles.insert(&c.id)
                || !sources.insert(&c.source_id)
                || !optional(&c.accession, 256)
                || c.description.len() > 4096
                || c.description.chars().any(|c| c.is_control() && c != '\t')
                || c.length == 0
                || c.length > MAX_BLOCKS as u64 * BLOCK_BASES as u64
                || !hex(&c.sequence_hash)
                || c.blocks.is_empty()
                || c.blocks.len() > MAX_BLOCKS
            {
                return Err(bad(
                    "invalid contig identity, length, description or block count",
                ));
            }
            let mut start = 0u64;
            for (bi, b) in c.blocks.iter().enumerate() {
                blocks = blocks.checked_add(1).ok_or(Error::TooLarge)?;
                if blocks > MAX_BLOCKS {
                    return Err(Error::TooLarge);
                }
                if b.index != bi as u64
                    || b.start != start
                    || b.count == 0
                    || b.count as usize > BLOCK_BASES
                    || b.offset != offset
                    || !hex(&b.encoded_hash)
                    || (bi + 1 < c.blocks.len() && b.count as usize != BLOCK_BASES)
                {
                    return Err(bad("block identity, order, size or file offset invalid"));
                }
                start = start
                    .checked_add(u64::from(b.count))
                    .ok_or(Error::TooLarge)?;
                offset = offset
                    .checked_add(u64::from(b.count).div_ceil(2))
                    .ok_or(Error::TooLarge)?;
            }
            if start != c.length {
                return Err(bad("block coverage differs from contig length"));
            }
            total = total.checked_add(c.length).ok_or(Error::TooLarge)?;
        }
        if total != self.total_bases
            || self.total_bases > MAX_BLOCKS as u64 * BLOCK_BASES as u64
            || self.content_revision != content_revision(&self.contigs)?
        {
            return Err(bad("canonical collection identity/count mismatch"));
        }
        let absent = missing(&self.build, &self.contigs)?;
        if self.missing_expected != absent
            || self.expected_complete
                != (!self.build.expected_components.is_empty() && absent.is_empty())
        {
            return Err(bad("expected-component completeness is inconsistent"));
        }
        encode(self)?;
        Ok(())
    }
    pub fn data_end(&self) -> u64 {
        self.contigs
            .last()
            .and_then(|c| c.blocks.last())
            .map_or(HEADER.len() as u64, |b| {
                b.offset.saturating_add(u64::from(b.count).div_ceil(2))
            })
    }
    pub fn manifest_hash(&self) -> Result<String, Error> {
        self.validate()?;
        Ok(blake3::hash(&encode(self)?).to_hex().to_string())
    }
}
pub fn encode_manifest(m: &GenomeManifest) -> Result<Vec<u8>, Error> {
    m.validate()?;
    encode(m)
}
impl FinishedReference {
    pub fn validate(&self) -> Result<(), Error> {
        self.manifest.validate()?;
        let bytes = encode(&self.manifest)?;
        if self.manifest_hash != blake3::hash(&bytes).to_hex().as_str()
            || self.manifest_offset != self.manifest.data_end()
            || self.file_length
                != self
                    .manifest_offset
                    .checked_add(bytes.len() as u64)
                    .and_then(|x| x.checked_add(TAIL_BYTES as u64))
                    .ok_or(Error::TooLarge)?
        {
            return Err(bad("finished layout/hash/file length mismatch"));
        }
        Ok(())
    }
}
pub fn encode_tail(f: &FinishedReference) -> Result<[u8; TAIL_BYTES], Error> {
    f.validate()?;
    let mut tail = [0u8; TAIL_BYTES];
    let len = f.file_length - f.manifest_offset - TAIL_BYTES as u64;
    tail[..8].copy_from_slice(&len.to_le_bytes());
    for (i, pair) in f.manifest_hash.as_bytes().chunks_exact(2).enumerate() {
        let digit = |c: u8| if c <= b'9' { c - b'0' } else { c - b'a' + 10 };
        tail[8 + i] = (digit(pair[0]) << 4) | digit(pair[1]);
    }
    tail[40..].copy_from_slice(&TAIL_MAGIC);
    Ok(tail)
}
pub fn parse_tail(header: &[u8], tail: &[u8], file_length: u64) -> Result<TailInfo, Error> {
    if header != HEADER
        || tail.len() != TAIL_BYTES
        || tail[40..] != TAIL_MAGIC
        || file_length < HEADER.len() as u64 + TAIL_BYTES as u64
    {
        return Err(bad("unsupported reference header/footer"));
    }
    let len = u64::from_le_bytes(tail[..8].try_into().unwrap());
    if len == 0 || len > MAX_MANIFEST_BYTES as u64 {
        return Err(Error::TooLarge);
    }
    let offset = file_length
        .checked_sub(TAIL_BYTES as u64)
        .and_then(|n| n.checked_sub(len))
        .filter(|n| *n >= HEADER.len() as u64)
        .ok_or_else(|| bad("manifest range outside file"))?;
    let hash = tail[8..40]
        .iter()
        .map(|b| format!("{b:02x}"))
        .collect::<String>();
    Ok(TailInfo {
        manifest_length: len,
        manifest_hash: hash,
        manifest_offset: offset,
        file_length,
    })
}
/// Decode only after a no-allocation typed Postcard shape preflight.
pub fn decode_manifest(bytes: &[u8], tail: &TailInfo) -> Result<FinishedReference, Error> {
    if bytes.len() > MAX_MANIFEST_BYTES
        || bytes.len() as u64 != tail.manifest_length
        || !hex(&tail.manifest_hash)
        || blake3::hash(bytes).to_hex().as_str() != tail.manifest_hash
    {
        return Err(bad("manifest byte length/hash mismatch"));
    }
    preflight_manifest(bytes)?;
    let (manifest, rest): (GenomeManifest, &[u8]) =
        postcard::take_from_bytes(bytes).map_err(|_| bad("malformed manifest encoding"))?;
    if !rest.is_empty() || encode(&manifest)? != bytes {
        return Err(bad("trailing or noncanonical manifest bytes"));
    }
    let finished = FinishedReference {
        manifest,
        manifest_hash: tail.manifest_hash.clone(),
        manifest_offset: tail.manifest_offset,
        file_length: tail.file_length,
    };
    finished.validate()?;
    Ok(finished)
}
pub fn encode_block(dna: &[u8]) -> Result<Vec<u8>, Error> {
    if dna.is_empty() || dna.len() > BLOCK_BASES {
        return Err(Error::TooLarge);
    }
    let mut out = vec![0u8; dna.len().div_ceil(2)];
    for (i, b) in dna.iter().enumerate() {
        if !b"ACGTRYSWKMBDHVN".contains(b) {
            return Err(bad("block contains invalid canonical DNA call"));
        }
        let mask = analysis::base_mask(*b);
        if i % 2 == 0 {
            out[i / 2] = mask << 4
        } else {
            out[i / 2] |= mask
        }
    }
    Ok(out)
}
fn decode_block(encoded: &[u8], count: u32) -> Result<Vec<u8>, Error> {
    let count = count as usize;
    if count == 0
        || count > BLOCK_BASES
        || encoded.len() != count.div_ceil(2)
        || count % 2 == 1 && encoded.last().is_some_and(|b| b & 15 != 0)
    {
        return Err(bad("packed block size/padding mismatch"));
    }
    let alphabet = b"?ACMGRSVTWYHKDBN";
    let mut dna = Vec::with_capacity(count);
    for i in 0..count {
        let mask = if i % 2 == 0 {
            encoded[i / 2] >> 4
        } else {
            encoded[i / 2] & 15
        };
        if mask == 0 {
            return Err(bad("zero IUPAC mask in valid block positions"));
        }
        dna.push(alphabet[mask as usize]);
    }
    Ok(dna)
}
struct Shape<'a> {
    bytes: &'a [u8],
    at: usize,
}
impl Shape<'_> {
    fn num(&mut self) -> Result<u64, Error> {
        let mut value = 0u64;
        for shift in (0..70).step_by(7) {
            let b = *self
                .bytes
                .get(self.at)
                .ok_or_else(|| bad("truncated manifest integer"))?;
            self.at += 1;
            if shift == 63 && b > 1 {
                return Err(bad("manifest integer overflow"));
            }
            value |= u64::from(b & 127) << shift;
            if b & 128 == 0 {
                if shift > 0 && b == 0 {
                    return Err(bad("noncanonical manifest integer"));
                }
                return Ok(value);
            }
        }
        Err(bad("manifest integer overflow"))
    }
    fn count(&mut self, max: usize) -> Result<usize, Error> {
        usize::try_from(self.num()?)
            .ok()
            .filter(|n| *n <= max)
            .ok_or(Error::TooLarge)
    }
    fn string(&mut self, max: usize) -> Result<(), Error> {
        let n = self.count(max)?;
        self.at = self
            .at
            .checked_add(n)
            .filter(|n| *n <= self.bytes.len())
            .ok_or_else(|| bad("truncated manifest string"))?;
        Ok(())
    }
    fn option(&mut self, max: usize) -> Result<(), Error> {
        match self.num()? {
            0 => Ok(()),
            1 => self.string(max),
            _ => Err(bad("invalid manifest option")),
        }
    }
    fn flag(&mut self) -> Result<(), Error> {
        if self.num()? > 1 {
            return Err(bad("invalid manifest enum/bool"));
        }
        Ok(())
    }
}
fn preflight_manifest(bytes: &[u8]) -> Result<(), Error> {
    if bytes.len() > MAX_MANIFEST_BYTES {
        return Err(Error::TooLarge);
    }
    let mut s = Shape { bytes, at: 0 };
    if s.num()? != 1 {
        return Err(bad("unsupported manifest schema"));
    }
    s.string(4096)?;
    s.option(256)?;
    s.option(256)?;
    for _ in 0..s.count(128)? {
        s.string(4096)?;
        s.option(512)?;
        s.option(128)?;
        s.option(64)?;
    }
    for _ in 0..s.count(MAX_CONTIGS)? {
        s.string(256)?;
        s.string(256)?;
        s.option(256)?;
        s.num()?;
        s.option(64)?;
        s.flag()?;
    }
    let mut total_blocks = 0usize;
    for _ in 0..s.count(MAX_CONTIGS)? {
        s.string(256)?;
        s.string(256)?;
        s.option(256)?;
        s.string(4096)?;
        s.num()?;
        s.flag()?;
        s.string(64)?;
        let count = s.count(MAX_BLOCKS - total_blocks)?;
        total_blocks += count;
        for _ in 0..count {
            s.num()?;
            s.num()?;
            if s.num()? > u64::from(u32::MAX) {
                return Err(Error::TooLarge);
            }
            s.num()?;
            s.string(64)?;
        }
    }
    s.string(64)?;
    s.string(64)?;
    s.num()?;
    s.num()?;
    for _ in 0..s.count(MAX_CONTIGS)? {
        s.string(256)?;
    }
    s.flag()?;
    if s.at != bytes.len() {
        return Err(bad("trailing manifest data"));
    }
    Ok(())
}
struct Current {
    record: ContigRecord,
    dna: Vec<u8>,
    hasher: blake3::Hasher,
}
pub struct FastaBuilder {
    info: BuildInfo,
    options: IngestOptions,
    contigs: Vec<ContigRecord>,
    current: Option<Current>,
    header: Vec<u8>,
    in_header: bool,
    line_start: bool,
    failed: bool,
    source_hasher: blake3::Hasher,
    source_bytes: u64,
    total_bases: u64,
    offset: u64,
    block_count: usize,
    manifest_budget: usize,
    source_ids: BTreeSet<String>,
    role_ids: BTreeSet<String>,
}
impl FastaBuilder {
    pub fn begin(info: BuildInfo, options: IngestOptions) -> Result<Self, Error> {
        info.validate()?;
        options.validate()?;
        let missing_budget = info
            .expected_components
            .iter()
            .map(|c| c.id.len() + 10)
            .sum::<usize>();
        let manifest_budget = encode(&info)?
            .len()
            .checked_add(missing_budget)
            .and_then(|n| n.checked_add(4096))
            .filter(|n| *n <= MAX_MANIFEST_BYTES)
            .ok_or(Error::TooLarge)?;
        Ok(Self {
            info,
            options,
            contigs: vec![],
            current: None,
            header: vec![],
            in_header: false,
            line_start: true,
            failed: false,
            source_hasher: blake3::Hasher::new(),
            source_bytes: 0,
            total_bases: 0,
            offset: HEADER.len() as u64,
            block_count: 0,
            manifest_budget,
            source_ids: BTreeSet::new(),
            role_ids: BTreeSet::new(),
        })
    }
    fn open_header(&mut self) -> Result<(), Error> {
        let raw = std::str::from_utf8(&self.header)
            .map_err(|_| bad("FASTA header is not UTF-8"))?
            .trim_end_matches('\r');
        if raw.is_empty() || raw.chars().any(|c| c.is_control() && c != '\t') {
            return Err(bad("FASTA header empty or contains controls"));
        }
        let split = raw.find(char::is_whitespace).unwrap_or(raw.len());
        let source_id = &raw[..split];
        if !id(source_id) || !self.source_ids.insert(source_id.to_owned()) {
            return Err(bad("FASTA source IDs must be unique bounded first tokens"));
        }
        if self.contigs.len() >= self.options.max_contigs {
            return Err(Error::TooLarge);
        }
        let expected = self
            .info
            .expected_components
            .iter()
            .find(|e| e.source_id == source_id);
        let role = expected.map_or(source_id, |e| e.id.as_str());
        if !self.role_ids.insert(role.to_owned()) {
            return Err(bad("contig role collision"));
        }
        let record = ContigRecord {
            id: role.into(),
            source_id: source_id.into(),
            accession: expected.and_then(|e| e.accession.clone()),
            description: raw[split..].trim().to_owned(),
            length: 0,
            topology: expected.map_or(Topology::Linear, |e| e.topology),
            sequence_hash: "0".repeat(64),
            blocks: vec![],
        };
        self.manifest_budget = self
            .manifest_budget
            .checked_add(encode(&record)?.len() + 16)
            .filter(|n| *n <= MAX_MANIFEST_BYTES)
            .ok_or(Error::TooLarge)?;
        self.current = Some(Current {
            record,
            dna: Vec::with_capacity(BLOCK_BASES),
            hasher: blake3::Hasher::new(),
        });
        self.header.clear();
        Ok(())
    }
    fn flush(&mut self, events: &mut Vec<EncodedBlock>) -> Result<(), Error> {
        let Some(current) = &mut self.current else {
            return Ok(());
        };
        if current.dna.is_empty() {
            return Ok(());
        }
        if self.block_count >= self.options.max_blocks {
            return Err(Error::TooLarge);
        }
        let encoded = encode_block(&current.dna)?;
        let descriptor = BlockDescriptor {
            index: current.record.blocks.len() as u64,
            start: current.record.length - current.dna.len() as u64,
            count: current.dna.len() as u32,
            offset: self.offset,
            encoded_hash: blake3::hash(&encoded).to_hex().to_string(),
        };
        self.manifest_budget = self
            .manifest_budget
            .checked_add(encode(&descriptor)?.len() + 2)
            .filter(|n| *n <= MAX_MANIFEST_BYTES)
            .ok_or(Error::TooLarge)?;
        current.hasher.update(&current.dna);
        current.dna.clear();
        self.offset = self
            .offset
            .checked_add(encoded.len() as u64)
            .ok_or(Error::TooLarge)?;
        self.block_count += 1;
        current.record.blocks.push(descriptor.clone());
        events.push(EncodedBlock {
            contig_id: current.record.id.clone(),
            descriptor,
            encoded,
        });
        Ok(())
    }
    fn close(&mut self, events: &mut Vec<EncodedBlock>) -> Result<(), Error> {
        self.flush(events)?;
        if let Some(mut c) = self.current.take() {
            if c.record.length == 0 {
                return Err(bad("empty FASTA contig"));
            }
            c.record.sequence_hash = c.hasher.finalize().to_hex().to_string();
            self.contigs.push(c.record);
        }
        Ok(())
    }
    pub fn push(&mut self, bytes: &[u8]) -> Result<Vec<EncodedBlock>, Error> {
        if self.failed {
            return Err(bad("ingestion session is poisoned; discard staging data"));
        }
        let result = self.push_inner(bytes);
        if result.is_err() {
            self.failed = true
        }
        result
    }
    fn push_inner(&mut self, bytes: &[u8]) -> Result<Vec<EncodedBlock>, Error> {
        if bytes.is_empty() || bytes.len() > self.options.max_push_bytes {
            return Err(Error::TooLarge);
        }
        self.source_bytes = self
            .source_bytes
            .checked_add(bytes.len() as u64)
            .filter(|n| *n <= self.options.max_source_bytes)
            .ok_or(Error::TooLarge)?;
        self.source_hasher.update(bytes);
        let mut events = vec![];
        for b in bytes {
            if self.in_header {
                if *b == b'\n' {
                    self.open_header()?;
                    self.in_header = false;
                    self.line_start = true;
                } else {
                    if self.header.len() >= 4096 {
                        return Err(Error::TooLarge);
                    }
                    self.header.push(*b);
                    self.line_start = false;
                }
                continue;
            }
            if *b == b'>' && self.line_start {
                self.close(&mut events)?;
                self.in_header = true;
                self.header.clear();
                self.line_start = false;
                continue;
            }
            if matches!(*b, b'\n' | b'\r') {
                self.line_start = true;
                continue;
            }
            if matches!(*b, b' ' | b'\t') {
                self.line_start = false;
                continue;
            }
            let canonical = b.to_ascii_uppercase();
            if !b"ACGTRYSWKMBDHVN".contains(&canonical) {
                return Err(bad(
                    "invalid FASTA DNA call; gaps, digits and comments are not silently converted",
                ));
            }
            let current = self
                .current
                .as_mut()
                .ok_or_else(|| bad("FASTA sequence before first header"))?;
            current.record.length = current
                .record
                .length
                .checked_add(1)
                .filter(|n| *n <= self.options.max_contig_bases)
                .ok_or(Error::TooLarge)?;
            self.total_bases = self
                .total_bases
                .checked_add(1)
                .filter(|n| *n <= self.options.max_total_bases)
                .ok_or(Error::TooLarge)?;
            current.dna.push(canonical);
            self.line_start = false;
            if current.dna.len() == BLOCK_BASES {
                self.flush(&mut events)?;
            }
        }
        Ok(events)
    }
    pub fn finish(mut self) -> Result<(Vec<EncodedBlock>, FinishedReference), Error> {
        if self.failed {
            return Err(bad("poisoned ingestion cannot publish a reference"));
        }
        let mut events = vec![];
        if self.in_header {
            self.open_header()?;
        }
        self.close(&mut events)?;
        if self.contigs.is_empty() {
            return Err(bad("no FASTA contigs"));
        }
        let absent = missing(&self.info, &self.contigs)?;
        let manifest = GenomeManifest {
            schema_version: 1,
            expected_complete: !self.info.expected_components.is_empty() && absent.is_empty(),
            missing_expected: absent,
            content_revision: content_revision(&self.contigs)?,
            build: self.info,
            contigs: self.contigs,
            raw_source_hash: self.source_hasher.finalize().to_hex().to_string(),
            source_bytes: self.source_bytes,
            total_bases: self.total_bases,
        };
        manifest.validate()?;
        let bytes = encode(&manifest)?;
        let finished = FinishedReference {
            manifest,
            manifest_hash: blake3::hash(&bytes).to_hex().to_string(),
            manifest_offset: self.offset,
            file_length: self
                .offset
                .checked_add(bytes.len() as u64)
                .and_then(|n| n.checked_add(TAIL_BYTES as u64))
                .ok_or(Error::TooLarge)?,
        };
        finished.validate()?;
        Ok((events, finished))
    }
}
/// A validated layout context caches metadata/hash once; it does not attest DNA.
pub struct BlockVerifier {
    manifest: GenomeManifest,
    manifest_hash: String,
}
pub struct VerifiedBlock {
    contig_id: String,
    descriptor: BlockDescriptor,
    dna: Vec<u8>,
    content_revision: String,
    manifest_hash: String,
}
impl BlockVerifier {
    pub fn begin(manifest: &GenomeManifest) -> Result<Self, Error> {
        manifest.validate()?;
        Ok(Self {
            manifest: manifest.clone(),
            manifest_hash: blake3::hash(&encode(manifest)?).to_hex().to_string(),
        })
    }
    pub fn verify(&self, role: &str, index: u64, encoded: &[u8]) -> Result<VerifiedBlock, Error> {
        let contig = self
            .manifest
            .contigs
            .iter()
            .find(|c| c.id == role)
            .ok_or_else(|| bad("unknown contig role"))?;
        let descriptor = usize::try_from(index)
            .ok()
            .and_then(|i| contig.blocks.get(i))
            .ok_or_else(|| bad("unknown contig block"))?;
        if encoded.len() != u64::from(descriptor.count).div_ceil(2) as usize
            || blake3::hash(encoded).to_hex().as_str() != descriptor.encoded_hash
        {
            return Err(bad("encoded block length/hash mismatch"));
        }
        let dna = decode_block(encoded, descriptor.count)?;
        Ok(VerifiedBlock {
            contig_id: role.into(),
            descriptor: descriptor.clone(),
            dna,
            content_revision: self.manifest.content_revision.clone(),
            manifest_hash: self.manifest_hash.clone(),
        })
    }
}
impl VerifiedBlock {
    pub fn verify(
        m: &GenomeManifest,
        role: &str,
        index: u64,
        encoded: &[u8],
    ) -> Result<Self, Error> {
        BlockVerifier::begin(m)?.verify(role, index, encoded)
    }
    pub fn dna(&self) -> &[u8] {
        &self.dna
    }
    pub fn descriptor(&self) -> &BlockDescriptor {
        &self.descriptor
    }
    pub fn contig_id(&self) -> &str {
        &self.contig_id
    }
    pub fn content_revision(&self) -> &str {
        &self.content_revision
    }
    pub fn manifest_hash(&self) -> &str {
        &self.manifest_hash
    }
}
struct ExpectedSequence {
    id: String,
    length: u64,
    sequence_hash: String,
    blocks: usize,
}
pub struct SemanticValidator {
    expected: Vec<ExpectedSequence>,
    manifest_hash: String,
    content_revision: String,
    contig: usize,
    next_block: u64,
    contig_bases: u64,
    total_bases: u64,
    hasher: blake3::Hasher,
    failed: bool,
}
impl SemanticValidator {
    pub fn begin(m: &GenomeManifest) -> Result<Self, Error> {
        m.validate()?;
        Ok(Self {
            expected: m
                .contigs
                .iter()
                .map(|c| ExpectedSequence {
                    id: c.id.clone(),
                    length: c.length,
                    sequence_hash: c.sequence_hash.clone(),
                    blocks: c.blocks.len(),
                })
                .collect(),
            manifest_hash: blake3::hash(&encode(m)?).to_hex().to_string(),
            content_revision: m.content_revision.clone(),
            contig: 0,
            next_block: 0,
            contig_bases: 0,
            total_bases: 0,
            hasher: blake3::Hasher::new(),
            failed: false,
        })
    }
    pub fn push(&mut self, block: &VerifiedBlock) -> Result<(), Error> {
        if self.failed {
            return Err(bad("semantic validator is poisoned"));
        }
        let result = self.push_inner(block);
        if result.is_err() {
            self.failed = true
        }
        result
    }
    fn push_inner(&mut self, b: &VerifiedBlock) -> Result<(), Error> {
        let expected = self
            .expected
            .get(self.contig)
            .ok_or_else(|| bad("extra block after semantic EOF"))?;
        if b.manifest_hash != self.manifest_hash
            || b.content_revision != self.content_revision
            || b.contig_id != expected.id
            || b.descriptor.index != self.next_block
            || b.descriptor.start != self.contig_bases
        {
            return Err(bad(
                "semantic blocks must match pinned layout and strict canonical order",
            ));
        }
        self.hasher.update(&b.dna);
        self.contig_bases = self
            .contig_bases
            .checked_add(b.dna.len() as u64)
            .ok_or(Error::TooLarge)?;
        self.total_bases = self
            .total_bases
            .checked_add(b.dna.len() as u64)
            .ok_or(Error::TooLarge)?;
        self.next_block += 1;
        if self.next_block == expected.blocks as u64 {
            if self.contig_bases != expected.length
                || self.hasher.finalize().to_hex().as_str() != expected.sequence_hash
            {
                return Err(bad(
                    "canonical contig DNA differs from manifest sequence hash",
                ));
            }
            self.contig += 1;
            self.next_block = 0;
            self.contig_bases = 0;
            self.hasher = blake3::Hasher::new();
        }
        Ok(())
    }
    pub fn finish(self) -> Result<SemanticReceipt, Error> {
        if self.failed || self.contig != self.expected.len() || self.next_block != 0 {
            return Err(bad(
                "full semantic verification has unverified or mismatched blocks",
            ));
        }
        Ok(SemanticReceipt {
            manifest_hash: self.manifest_hash,
            content_revision: self.content_revision,
            contigs_verified: self.contig,
            total_bases: self.total_bases,
        })
    }
}
impl GenomicInterval {
    pub fn validate_for(&self, m: &GenomeManifest) -> Result<(), Error> {
        let c = m
            .contigs
            .iter()
            .find(|c| c.id == self.contig_id)
            .ok_or_else(|| bad("unknown global interval contig role"))?;
        if self.start >= self.end || self.end > c.length {
            return Err(bad("invalid u64 half-open contig interval"));
        }
        Ok(())
    }
}
impl GlobalPosition {
    pub fn validate_for(&self, m: &GenomeManifest) -> Result<(), Error> {
        if self.content_revision != m.content_revision {
            return Err(bad("stale global position reference identity"));
        }
        let c = m
            .contigs
            .iter()
            .find(|c| c.id == self.contig_id)
            .ok_or_else(|| bad("unknown global position contig"))?;
        if self.position >= c.length {
            return Err(bad("global base position outside contig"));
        }
        Ok(())
    }
}
impl ReferenceLink {
    pub fn validate(&self) -> Result<(), Error> {
        if self.schema_version != 1
            || !id(&self.role)
            || !hex(&self.manifest_hash)
            || !hex(&self.content_revision)
            || !text(&self.name, 4096)
            || !optional(&self.assembly_accession, 256)
            || !optional(&self.version, 256)
            || self.contig_count == 0
            || self.contig_count > MAX_CONTIGS
            || self.total_bases == 0
            || self.total_bases > MAX_BLOCKS as u64 * BLOCK_BASES as u64
            || self.file_length < HEADER.len() as u64 + TAIL_BYTES as u64
            || self.file_length
                > MAX_BLOCKS as u64 * (BLOCK_BASES as u64).div_ceil(2)
                    + MAX_MANIFEST_BYTES as u64
                    + TAIL_BYTES as u64
                    + HEADER.len() as u64
            || !optional(&self.cache_hint, 1024)
        {
            return Err(bad(
                "invalid reference link; it is descriptive identity, not data/search authority",
            ));
        }
        Ok(())
    }
}
pub fn reference_link(
    f: &FinishedReference,
    role: &str,
    cache_hint: Option<String>,
) -> Result<ReferenceLink, Error> {
    f.validate()?;
    let link = ReferenceLink {
        schema_version: 1,
        role: role.into(),
        manifest_hash: f.manifest_hash.clone(),
        content_revision: f.manifest.content_revision.clone(),
        name: f.manifest.build.name.clone(),
        assembly_accession: f.manifest.build.assembly_accession.clone(),
        version: f.manifest.build.version.clone(),
        contig_count: f.manifest.contigs.len(),
        total_bases: f.manifest.total_bases,
        file_length: f.file_length,
        cache_hint,
    };
    link.validate()?;
    Ok(link)
}
pub fn make_attachment(
    link: &ReferenceLink,
) -> Result<(String, crate::documents::workspace::WorkspaceAttachment), Error> {
    link.validate()?;
    let bytes = serde_json::to_vec(link).map_err(|_| Error::TooLarge)?;
    if bytes.len() > 32 * 1024 {
        return Err(Error::TooLarge);
    }
    let hash = blake3::hash(&bytes).to_hex().to_string();
    Ok((
        hash.clone(),
        crate::documents::workspace::WorkspaceAttachment {
            kind: WORKSPACE_ATTACHMENT_KIND.into(),
            name: format!("genome-reference-{hash}.json"),
            media_type: "application/json".into(),
            bytes,
        },
    ))
}
pub fn load_link(bytes: &[u8]) -> Result<ReferenceLink, Error> {
    if bytes.len() > 32 * 1024 {
        return Err(Error::TooLarge);
    }
    let link: ReferenceLink =
        serde_json::from_slice(bytes).map_err(|_| bad("invalid strict reference-link JSON"))?;
    link.validate()?;
    Ok(link)
}
