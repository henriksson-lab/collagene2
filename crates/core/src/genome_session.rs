//! Private live genome parser/verification/search state. Dataset bytes live with the host.
//! Serializable progress is evidence only and cannot be used to restore a live cursor.
use crate::{
    genome::{self, *},
    genome_crispr::{self, *},
};
use collagene_model::Error;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
macro_rules! dto{($($x:tt)*)=>{#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($x)*};}
dto! {pub struct GenomeRef{pub id:String,pub manifest_hash:String,pub content_revision:String}}
dto! {pub struct ImportStarted{pub binding:String,pub import_id:String,pub reference:Option<FinishedReference>,pub next_block:Option<BlockRequest>}}
dto! {pub struct ImportProgress{pub import_id:String,pub received:u64,pub blocks:Vec<EncodedBlock>,pub next_block:Option<BlockRequest>}}
dto! {pub struct FastaFinished{pub import_id:String,pub blocks:Vec<EncodedBlock>,pub reference:FinishedReference}}
dto! {pub struct RegisteredGenome{pub binding:String,pub genome:GenomeRef,pub manifest:GenomeManifest,pub semantic_verified:bool,pub receipt:Option<SemanticReceipt>}}
dto! {pub struct SearchState{pub binding:String,pub job_id:String,pub genome:GenomeRef,pub progress:GenomeSearchProgress}}
dto! {pub struct SearchPage{pub job_id:String,pub genome:GenomeRef,pub progress:GenomeSearchProgress,pub offset:usize,pub next_offset:usize,pub done:bool,pub hits:Vec<GenomeOffTargetHit>}}
dto! {pub struct IntervalGuides{pub genome:GenomeRef,pub report:GenomeGuideReport}}
dto! {pub struct IntervalRead{pub genome:GenomeRef,pub report:GenomeIntervalRead}}
dto! {pub struct Disposed{pub id:String,pub removed:bool}}
use std::sync::atomic::{AtomicU64, Ordering};
static NEXT_SESSION: AtomicU64 = AtomicU64::new(1);
const MAX_IMPORTS: usize = 2;
const MAX_DATASETS: usize = 4;
const MAX_SEARCHES: usize = 4;
fn bad(s: &str) -> Error {
    Error::Format(format!("genome session: {s}"))
}
fn token(s: &str) -> bool {
    !s.is_empty() && s.len() <= 256 && !s.chars().any(char::is_control)
}
struct Bounded(Vec<u8>);
impl postcard::ser_flavors::Flavor for Bounded {
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
pub fn fasta_revision(
    build: &BuildInfo,
    options: &IngestOptions,
    expected_source_bytes: u64,
) -> Result<String, Error> {
    build.validate()?;
    options.validate()?;
    if expected_source_bytes == 0 || expected_source_bytes > options.max_source_bytes {
        return Err(bad("declared source byte length exceeds import bounds"));
    }
    let b = postcard::serialize_with_flavor(
        &("genome-fasta-v1", build, options, expected_source_bytes),
        Bounded(vec![]),
    )
    .map_err(|_| Error::TooLarge)?;
    Ok(blake3::hash(&b).to_hex().to_string())
}
pub fn native_revision(
    header: &[u8],
    tail: &[u8],
    file_length: u64,
    manifest: &[u8],
) -> Result<String, Error> {
    if manifest.len() > MAX_MANIFEST_BYTES {
        return Err(Error::TooLarge);
    }
    let info = genome::parse_tail(header, tail, file_length)?;
    if manifest.len() as u64 != info.manifest_length
        || blake3::hash(manifest).to_hex().as_str() != info.manifest_hash
    {
        return Err(bad("native manifest bytes do not match footer"));
    }
    Ok(info.manifest_hash)
}
struct FastaUpload {
    builder: FastaBuilder,
    received: u64,
    max_push: usize,
    expected_source_bytes: u64,
}
struct NativeValidation {
    reference: FinishedReference,
    verifier: BlockVerifier,
    semantic: SemanticValidator,
    ci: usize,
    bi: usize,
    received: u64,
}
impl NativeValidation {
    fn next(&self) -> Option<BlockRequest> {
        self.reference
            .manifest
            .contigs
            .get(self.ci)
            .map(|c| BlockRequest {
                contig_id: c.id.clone(),
                index: self.bi as u64,
            })
    }
}
struct Dataset {
    manifest: GenomeManifest,
    reference: GenomeRef,
    receipt: Option<SemanticReceipt>,
}
struct Search {
    genome: GenomeRef,
    state: GenomeSearchSession,
}
pub struct GenomeSession {
    epoch: u64,
    serial: u64,
    fasta: BTreeMap<String, FastaUpload>,
    native: BTreeMap<String, NativeValidation>,
    datasets: BTreeMap<String, Dataset>,
    searches: BTreeMap<String, Search>,
}
impl Default for GenomeSession {
    fn default() -> Self {
        Self {
            epoch: NEXT_SESSION.fetch_add(1, Ordering::Relaxed),
            serial: 0,
            fasta: BTreeMap::new(),
            native: BTreeMap::new(),
            datasets: BTreeMap::new(),
            searches: BTreeMap::new(),
        }
    }
}
impl GenomeSession {
    fn id(&mut self, kind: &str) -> Result<String, Error> {
        self.serial = self.serial.checked_add(1).ok_or(Error::TooLarge)?;
        Ok(format!("g{}-{kind}-{}", self.epoch, self.serial))
    }
    fn dataset(&self, r: &GenomeRef) -> Result<&Dataset, Error> {
        let d = self
            .datasets
            .get(&r.id)
            .ok_or_else(|| bad("dataset handle is not live in this worker session"))?;
        if &d.reference != r {
            return Err(bad("dataset manifest/content binding changed"));
        }
        Ok(d)
    }
    fn register(
        &mut self,
        manifest: GenomeManifest,
        receipt: Option<SemanticReceipt>,
        binding: String,
    ) -> Result<RegisteredGenome, Error> {
        manifest.validate()?;
        let hash = manifest.manifest_hash()?;
        if let Some(proof) = &receipt {
            if proof.manifest_hash != hash
                || proof.content_revision != manifest.content_revision
                || proof.total_bases != manifest.total_bases
                || proof.contigs_verified != manifest.contigs.len()
            {
                return Err(bad("semantic receipt binding mismatch"));
            }
        }
        if receipt.is_some() {
            if let Some(d) = self
                .datasets
                .values_mut()
                .find(|d| d.reference.manifest_hash == hash)
            {
                d.receipt = receipt.clone();
                return Ok(RegisteredGenome {
                    binding,
                    genome: d.reference.clone(),
                    manifest: d.manifest.clone(),
                    semantic_verified: true,
                    receipt,
                });
            }
        }
        if self.datasets.len() >= MAX_DATASETS {
            return Err(bad("at most four live genome manifests"));
        }
        let reference = GenomeRef {
            id: self.id("dataset")?,
            manifest_hash: hash,
            content_revision: manifest.content_revision.clone(),
        };
        self.datasets.insert(
            reference.id.clone(),
            Dataset {
                manifest: manifest.clone(),
                reference: reference.clone(),
                receipt: receipt.clone(),
            },
        );
        Ok(RegisteredGenome {
            binding,
            genome: reference,
            manifest,
            semantic_verified: receipt.is_some(),
            receipt,
        })
    }
    pub fn process(
        &mut self,
        job: crate::protocol::schema::AnalysisJob,
        revision: &str,
    ) -> Result<crate::protocol::schema::AnalysisResult, Error> {
        use crate::protocol::schema::{AnalysisJob as J, AnalysisResult as R};
        match job {
            J::BeginGenomeFasta {
                build,
                options,
                expected_source_bytes,
            } => {
                if self.fasta.len() + self.native.len() >= MAX_IMPORTS {
                    return Err(bad("at most two live genome imports"));
                }
                if fasta_revision(&build, &options, expected_source_bytes)? != revision {
                    return Err(bad("FASTA options/build binding changed"));
                }
                let max_push = options.max_push_bytes;
                let builder = FastaBuilder::begin(build, options)?;
                let id = self.id("fasta")?;
                self.fasta.insert(
                    id.clone(),
                    FastaUpload {
                        builder,
                        received: 0,
                        max_push,
                        expected_source_bytes,
                    },
                );
                Ok(R::GenomeImportStarted(ImportStarted {
                    binding: revision.into(),
                    import_id: id,
                    reference: None,
                    next_block: None,
                }))
            }
            J::PushGenomeFasta {
                import_id,
                offset,
                bytes,
            } => {
                if revision != import_id {
                    return Err(bad("import binding mismatch"));
                }
                let current = self
                    .fasta
                    .get(&import_id)
                    .ok_or_else(|| bad("FASTA import is not live"))?;
                if offset != current.received {
                    return Err(bad("FASTA offset mismatch"));
                }
                if bytes.is_empty() || bytes.len() > current.max_push {
                    return Err(Error::TooLarge);
                }
                let received = offset
                    .checked_add(bytes.len() as u64)
                    .ok_or(Error::TooLarge)?;
                if received > current.expected_source_bytes {
                    return Err(bad("FASTA chunk exceeds declared source length"));
                }
                let mut upload = self.fasta.remove(&import_id).unwrap();
                let blocks = upload.builder.push(&bytes)?;
                upload.received = received;
                self.fasta.insert(import_id.clone(), upload);
                Ok(R::GenomeImportProgress(ImportProgress {
                    import_id,
                    received,
                    blocks,
                    next_block: None,
                }))
            }
            J::FinishGenomeFasta { import_id } => {
                if revision != import_id {
                    return Err(bad("import binding mismatch"));
                }
                let current = self
                    .fasta
                    .get(&import_id)
                    .ok_or_else(|| bad("FASTA import is not live"))?;
                if current.received != current.expected_source_bytes {
                    return Err(bad(
                        "FASTA source is truncated relative to declared byte length",
                    ));
                }
                let upload = self
                    .fasta
                    .remove(&import_id)
                    .ok_or_else(|| bad("FASTA import is not live"))?;
                let (blocks, reference) = upload.builder.finish()?;
                Ok(R::GenomeFastaFinished(FastaFinished {
                    import_id,
                    blocks,
                    reference,
                }))
            }
            J::BeginGenomeNative {
                header,
                tail,
                file_length,
                manifest,
            } => {
                if self.fasta.len() + self.native.len() >= MAX_IMPORTS {
                    return Err(bad("at most two live genome imports"));
                }
                if native_revision(&header, &tail, file_length, &manifest)? != revision {
                    return Err(bad("native manifest binding mismatch"));
                }
                let info = genome::parse_tail(&header, &tail, file_length)?;
                let reference = genome::decode_manifest(&manifest, &info)?;
                let verifier = BlockVerifier::begin(&reference.manifest)?;
                let semantic = SemanticValidator::begin(&reference.manifest)?;
                let id = self.id("native")?;
                let validation = NativeValidation {
                    reference: reference.clone(),
                    verifier,
                    semantic,
                    ci: 0,
                    bi: 0,
                    received: 0,
                };
                let next_block = validation.next();
                self.native.insert(id.clone(), validation);
                Ok(R::GenomeImportStarted(ImportStarted {
                    binding: revision.into(),
                    import_id: id,
                    reference: Some(reference),
                    next_block,
                }))
            }
            J::PushGenomeNativeBlock { import_id, block } => {
                if revision != import_id {
                    return Err(bad("import binding mismatch"));
                }
                if block.encoded.len() > BLOCK_BASES.div_ceil(2) {
                    return Err(Error::TooLarge);
                }
                let current = self
                    .native
                    .get(&import_id)
                    .ok_or_else(|| bad("native validator is not live"))?;
                let next = current
                    .next()
                    .ok_or_else(|| bad("all native blocks already verified"))?;
                if next.contig_id != block.contig_id || next.index != block.descriptor.index {
                    return Err(bad("native blocks must follow exact manifest order"));
                }
                let verified = current.verifier.verify(
                    &block.contig_id,
                    block.descriptor.index,
                    &block.encoded,
                )?;
                if verified.descriptor() != &block.descriptor {
                    return Err(bad("block descriptor differs from manifest"));
                }
                let mut validation = self.native.remove(&import_id).unwrap();
                validation.semantic.push(&verified)?;
                validation.received = validation
                    .received
                    .checked_add(u64::from(block.descriptor.count))
                    .ok_or(Error::TooLarge)?;
                validation.bi += 1;
                if validation.bi
                    == validation.reference.manifest.contigs[validation.ci]
                        .blocks
                        .len()
                {
                    validation.ci += 1;
                    validation.bi = 0;
                }
                let received = validation.received;
                let next_block = validation.next();
                self.native.insert(import_id.clone(), validation);
                Ok(R::GenomeImportProgress(ImportProgress {
                    import_id,
                    received,
                    blocks: vec![],
                    next_block,
                }))
            }
            J::FinishGenomeNative { import_id } => {
                if revision != import_id {
                    return Err(bad("import binding mismatch"));
                }
                let current = self
                    .native
                    .get(&import_id)
                    .ok_or_else(|| bad("native validator is not live"))?;
                if current.next().is_some() {
                    return Err(bad(
                        "native semantic validation has not reached every block",
                    ));
                }
                let hash = &current.reference.manifest_hash;
                if self.datasets.len() >= MAX_DATASETS
                    && !self
                        .datasets
                        .values()
                        .any(|d| &d.reference.manifest_hash == hash)
                {
                    return Err(bad("at most four live genome manifests"));
                }
                let validation = self.native.remove(&import_id).unwrap();
                let receipt = validation.semantic.finish()?;
                self.register(
                    validation.reference.manifest,
                    Some(receipt),
                    revision.into(),
                )
                .map(R::GenomeRegistered)
            }
            J::AbortGenomeImport { import_id } => {
                if revision != import_id {
                    return Err(bad("import binding mismatch"));
                }
                let removed = self.fasta.remove(&import_id).is_some()
                    | self.native.remove(&import_id).is_some();
                Ok(R::GenomeDisposed(Disposed {
                    id: import_id,
                    removed,
                }))
            }
            J::RegisterGenome { manifest } => {
                if manifest.manifest_hash()? != revision {
                    return Err(bad("registration manifest binding mismatch"));
                }
                self.register(manifest, None, revision.into())
                    .map(R::GenomeRegistered)
            }
            J::GenomeIntervalGuides {
                genome,
                options,
                blocks,
            } => {
                if genome.manifest_hash != revision {
                    return Err(bad("guide dataset binding mismatch"));
                }
                let d = self.dataset(&genome)?;
                let report = genome_crispr::design_interval(&d.manifest, &options, &blocks)?;
                Ok(R::GenomeIntervalGuides(IntervalGuides { genome, report }))
            }
            J::ReadGenomeInterval {
                genome,
                options,
                blocks,
            } => {
                if genome.manifest_hash != revision {
                    return Err(bad("window dataset binding mismatch"));
                }
                let d = self.dataset(&genome)?;
                let report = genome_crispr::read_interval(&d.manifest, &options, &blocks)?;
                Ok(R::GenomeIntervalRead(IntervalRead { genome, report }))
            }
            J::StartGenomeSearch {
                genome,
                query,
                options,
            } => {
                if genome.manifest_hash != revision {
                    return Err(bad("search dataset binding mismatch"));
                }
                if self.searches.len() >= MAX_SEARCHES {
                    return Err(bad("at most four live genome searches"));
                }
                let d = self.dataset(&genome)?;
                let state = GenomeSearchSession::begin(&d.manifest, &query, &options)?;
                let progress = state.progress();
                let id = self.id("search")?;
                self.searches.insert(
                    id.clone(),
                    Search {
                        genome: genome.clone(),
                        state,
                    },
                );
                Ok(R::GenomeSearchState(SearchState {
                    binding: revision.into(),
                    job_id: id,
                    genome,
                    progress,
                }))
            }
            J::AdvanceGenomeSearch { job_id, blocks } => {
                if revision != job_id {
                    return Err(bad("search handle binding mismatch"));
                }
                let search = self
                    .searches
                    .get_mut(&job_id)
                    .ok_or_else(|| bad("search is not live"))?;
                let progress = search.state.advance(&blocks)?;
                Ok(R::GenomeSearchState(SearchState {
                    binding: revision.into(),
                    job_id,
                    genome: search.genome.clone(),
                    progress,
                }))
            }
            J::PageGenomeSearch {
                job_id,
                offset,
                limit,
            } => {
                if revision != job_id {
                    return Err(bad("search handle binding mismatch"));
                }
                let search = self
                    .searches
                    .get(&job_id)
                    .ok_or_else(|| bad("search is not live"))?;
                let hits = search.state.page(offset, limit)?;
                let progress = search.state.progress();
                let next_offset = offset.checked_add(hits.len()).ok_or(Error::TooLarge)?;
                let done = next_offset == progress.hits_retained;
                Ok(R::GenomeSearchPage(SearchPage {
                    job_id,
                    genome: search.genome.clone(),
                    progress,
                    offset,
                    next_offset,
                    done,
                    hits,
                }))
            }
            J::CancelGenomeJob { id } => {
                if revision != id {
                    return Err(bad("resource handle binding mismatch"));
                }
                let mut removed = self.searches.remove(&id).is_some();
                removed |= self.fasta.remove(&id).is_some();
                removed |= self.native.remove(&id).is_some();
                if self.datasets.remove(&id).is_some() {
                    self.searches.retain(|_, s| s.genome.id != id);
                    removed = true;
                }
                Ok(R::GenomeDisposed(Disposed { id, removed }))
            }
            _ => Err(bad("job requires another worker dispatcher")),
        }
    }
}
fn digest(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
impl GenomeRef {
    pub fn validate(&self) -> Result<(), Error> {
        if !token(&self.id) || !digest(&self.manifest_hash) || !digest(&self.content_revision) {
            return Err(bad("invalid dataset handle binding"));
        }
        Ok(())
    }
}
fn block_shape(b: &EncodedBlock) -> Result<(), Error> {
    let d = &b.descriptor;
    if !token(&b.contig_id)
        || d.count == 0
        || d.count as usize > BLOCK_BASES
        || d.start
            != d.index
                .checked_mul(BLOCK_BASES as u64)
                .ok_or(Error::TooLarge)?
        || b.encoded.len() != (d.count as usize).div_ceil(2)
        || !digest(&d.encoded_hash)
        || blake3::hash(&b.encoded).to_hex().as_str() != d.encoded_hash
    {
        return Err(bad("invalid bounded block transport"));
    }
    Ok(())
}
fn proof_shape(
    manifest: &GenomeManifest,
    reference: &GenomeRef,
    proof: Option<&SemanticReceipt>,
    verified: bool,
) -> Result<(), Error> {
    reference.validate()?;
    manifest.validate()?;
    if manifest.manifest_hash()? != reference.manifest_hash
        || manifest.content_revision != reference.content_revision
        || verified != proof.is_some()
    {
        return Err(bad("manifest/receipt transport binding mismatch"));
    }
    if let Some(p) = proof {
        if p.manifest_hash != reference.manifest_hash
            || p.content_revision != reference.content_revision
            || p.contigs_verified != manifest.contigs.len()
            || p.total_bases != manifest.total_bases
        {
            return Err(bad("semantic receipt transport mismatch"));
        }
    }
    Ok(())
}
/// Lightweight transport binding. This never restores private parser/search authority.
pub fn validate_result(
    result: &crate::protocol::schema::AnalysisResult,
    revision: &str,
) -> Result<(), Error> {
    use crate::protocol::schema::AnalysisResult as R;
    match result {
        R::GenomeImportStarted(r) => {
            if r.binding != revision || !token(&r.import_id) {
                return Err(bad("import receipt binding mismatch"));
            }
            if let Some(f) = &r.reference {
                f.validate()?;
                if f.manifest_hash != revision {
                    return Err(bad("native import source mismatch"));
                }
                if let Some(next) = &r.next_block {
                    let c = f
                        .manifest
                        .contigs
                        .first()
                        .ok_or_else(|| bad("missing first contig"))?;
                    if next.contig_id != c.id || next.index != 0 {
                        return Err(bad("native import initial cursor mismatch"));
                    }
                }
            } else if r.next_block.is_some() {
                return Err(bad("FASTA import has unexpected native cursor"));
            }
        }
        R::GenomeImportProgress(r) => {
            if r.import_id != revision || !token(&r.import_id) || r.blocks.len() > 10_016 {
                return Err(bad("import progress binding mismatch"));
            }
            for b in &r.blocks {
                block_shape(b)?;
            }
            if let Some(next) = &r.next_block {
                if !token(&next.contig_id) {
                    return Err(bad("invalid next native block"));
                }
            }
        }
        R::GenomeFastaFinished(r) => {
            if r.import_id != revision || !token(&r.import_id) || r.blocks.len() > 10_016 {
                return Err(bad("FASTA finish binding mismatch"));
            }
            r.reference.validate()?;
            for b in &r.blocks {
                block_shape(b)?;
            }
        }
        R::GenomeRegistered(r) => {
            if r.binding != revision {
                return Err(bad("registration response binding mismatch"));
            }
            proof_shape(
                &r.manifest,
                &r.genome,
                r.receipt.as_ref(),
                r.semantic_verified,
            )?;
        }
        R::GenomeIntervalGuides(r) => {
            r.genome.validate()?;
            r.report.validate()?;
            if r.genome.manifest_hash != revision
                || r.report.manifest_hash != r.genome.manifest_hash
                || r.report.content_revision != r.genome.content_revision
            {
                return Err(bad("interval guide response mismatch"));
            }
        }
        R::GenomeIntervalRead(r) => {
            r.genome.validate()?;
            r.report.validate()?;
            if r.genome.manifest_hash != revision
                || r.report.manifest_hash != r.genome.manifest_hash
                || r.report.content_revision != r.genome.content_revision
            {
                return Err(bad("interval read response mismatch"));
            }
        }
        R::GenomeSearchState(r) => {
            r.genome.validate()?;
            r.progress.validate()?;
            if r.binding != revision
                || !token(&r.job_id)
                || r.progress.manifest_hash != r.genome.manifest_hash
                || r.progress.content_revision != r.genome.content_revision
            {
                return Err(bad("search state response mismatch"));
            }
        }
        R::GenomeSearchPage(r) => {
            r.genome.validate()?;
            r.progress.validate()?;
            if r.job_id != revision
                || !token(&r.job_id)
                || r.hits.len() > 256
                || r.next_offset != r.offset.checked_add(r.hits.len()).ok_or(Error::TooLarge)?
                || r.next_offset > r.progress.hits_retained
                || r.done != (r.next_offset == r.progress.hits_retained)
                || r.progress.manifest_hash != r.genome.manifest_hash
                || r.progress.content_revision != r.genome.content_revision
            {
                return Err(bad("search page response mismatch"));
            }
            for h in &r.hits {
                if !digest(&h.id)
                    || !token(&h.contig_id)
                    || h.guide_positions.len() != r.progress.options.profile.guide_length
                    || h.pam_positions.len() != r.progress.options.profile.pam.len()
                    || h.reference_guide.len() != h.guide_positions.len()
                    || h.pam_sequence.len() != h.pam_positions.len()
                    || h.min_mismatches > h.max_mismatches
                    || h.max_mismatches > h.guide_positions.len()
                {
                    return Err(bad("invalid bounded search hit page"));
                }
            }
        }
        R::GenomeDisposed(r) => {
            if r.id != revision || !token(&r.id) {
                return Err(bad("disposal response binding mismatch"));
            }
        }
        _ => return Err(bad("not a genome session response")),
    }
    Ok(())
}
