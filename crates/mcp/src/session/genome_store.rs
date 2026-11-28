//! Opaque session-owned reference files. No caller-controlled filesystem paths.

use collagene_core::tracks::{
    crispr::BlockRequest,
    genome::{
        self, BlockVerifier, BuildInfo, EncodedBlock, FastaBuilder, FinishedReference,
        IngestOptions, SemanticValidator,
    },
};
use collagene_model::Error;
use serde::{Deserialize, Serialize};
use std::{
    collections::BTreeMap,
    fs::File,
    io::{Read, Seek, SeekFrom, Write},
};

const MAX_IMPORTS: usize = 2;
const MAX_DATASETS: usize = 4;
const MAX_MANAGED_BYTES: u64 = 64 * 1024 * 1024 * 1024;
const MAX_CHUNK: usize = 1024 * 1024;
use collagene_local_storage::{ManagedFile, ManagedSession, StorageRoot};
#[cfg(test)]
use std::{fs, path::PathBuf};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, schemars::JsonSchema)]
pub enum ImportFormat {
    Fasta,
    Cgref,
}
#[derive(Debug, Clone, Serialize)]
pub(crate) struct NativeDataset {
    pub id: String,
    pub reference: FinishedReference,
    pub semantic_verified: bool,
    pub input_bytes: u64,
}
#[derive(Debug, Clone, Serialize)]
pub(crate) struct NativeImportProgress {
    pub done: bool,
    pub import_id: String,
    pub source_bytes: u64,
    pub blocks_verified: u64,
    pub blocks_total: u64,
    pub dataset: Option<NativeDataset>,
}
fn bad(s: &str) -> Error {
    Error::Format(format!("native reference store: {s}"))
}
fn io(e: std::io::Error) -> Error {
    Error::Unsupported(format!("managed reference I/O: {e}"))
}
fn digest(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_hexdigit() && !b.is_ascii_uppercase())
}
type OwnedFile = ManagedFile;
struct Verification {
    reference: FinishedReference,
    verifier: BlockVerifier,
    semantic: Option<SemanticValidator>,
    contig: usize,
    block: usize,
    verified: u64,
    total: u64,
}
struct Import {
    format: ImportFormat,
    expected: u64,
    received: u64,
    expected_hash: Option<String>,
    input_hasher: blake3::Hasher,
    reserved: u64,
    file: OwnedFile,
    fasta: Option<FastaBuilder>,
    verification: Option<Verification>,
}
struct Dataset {
    info: NativeDataset,
    file: OwnedFile,
}
pub(crate) struct NativeGenomes {
    #[cfg(test)]
    directory: PathBuf,
    session: ManagedSession,
    namespace: String,
    serial: u64,
    imports: BTreeMap<String, Import>,
    datasets: BTreeMap<String, Dataset>,
    completed: BTreeMap<String, String>,
}
impl NativeGenomes {
    pub fn new() -> Result<Self, Error> {
        let root = StorageRoot::temporary().map_err(io)?;
        Self::in_root(&root).map_err(io)
    }
    pub fn in_root(root: &StorageRoot) -> std::io::Result<Self> {
        let session = root.session()?;
        Ok(Self {
            #[cfg(test)]
            directory: session.diagnostic_path(),
            namespace: session.namespace().to_owned(),
            session,
            serial: 0,
            imports: BTreeMap::new(),
            datasets: BTreeMap::new(),
            completed: BTreeMap::new(),
        })
    }
    fn usage(&self) -> u64 {
        self.imports.values().map(|s| s.reserved).sum::<u64>()
            + self
                .datasets
                .values()
                .map(|d| d.info.reference.file_length)
                .sum::<u64>()
    }
    pub fn begin(
        &mut self,
        format: ImportFormat,
        expected_bytes: u64,
        expected_hash: Option<String>,
        build: Option<BuildInfo>,
        options: Option<IngestOptions>,
    ) -> Result<String, Error> {
        if self.imports.len() >= MAX_IMPORTS {
            return Err(bad("at most two active reference imports"));
        }
        if expected_bytes == 0
            || expected_bytes > MAX_MANAGED_BYTES
            || expected_hash.as_ref().is_some_and(|h| !digest(h))
        {
            return Err(bad("invalid expected input length/hash"));
        }
        if format == ImportFormat::Cgref && (build.is_some() || options.is_some()) {
            return Err(bad(
                "native reference imports preserve their manifest; build/options overrides are forbidden",
            ));
        }
        let reserved = expected_bytes
            .checked_add(if format == ImportFormat::Fasta {
                genome::MAX_MANIFEST_BYTES as u64
                    + genome::HEADER.len() as u64
                    + genome::TAIL_BYTES as u64
            } else {
                0
            })
            .filter(|n| *n <= MAX_MANAGED_BYTES.saturating_sub(self.usage()))
            .ok_or(Error::TooLarge)?;
        let fasta = if format == ImportFormat::Fasta {
            let options = options.unwrap_or_default();
            options.validate()?;
            if expected_bytes > options.max_source_bytes {
                return Err(Error::TooLarge);
            }
            Some(FastaBuilder::begin(build.unwrap_or_default(), options)?)
        } else {
            None
        };
        self.serial = self.serial.checked_add(1).ok_or(Error::TooLarge)?;
        let id = format!("reference-import-{}-{}", self.namespace, self.serial);
        let mut owned = self.session.create_file().map_err(io)?;
        if format == ImportFormat::Fasta {
            owned.file().write_all(&genome::HEADER).map_err(io)?;
        }
        self.imports.insert(
            id.clone(),
            Import {
                format,
                expected: expected_bytes,
                received: 0,
                expected_hash,
                input_hasher: blake3::Hasher::new(),
                reserved,
                file: owned,
                fasta,
                verification: None,
            },
        );
        Ok(id)
    }
    pub fn append(&mut self, id: &str, offset: u64, bytes: &[u8]) -> Result<u64, Error> {
        let stage = self
            .imports
            .get(id)
            .ok_or_else(|| bad("unknown reference import"))?;
        if bytes.is_empty()
            || bytes.len() > MAX_CHUNK
            || offset != stage.received
            || stage.verification.is_some()
            || bytes.len() as u64 > stage.expected.saturating_sub(stage.received)
        {
            return Err(bad(
                "append needs exact next offset and 1..1MiB bytes within expected input",
            ));
        }
        let mut stage = self.imports.remove(id).unwrap();
        let result = (|| {
            if let Some(builder) = &mut stage.fasta {
                let events = builder.push(bytes)?;
                write_events(&mut stage.file, &events)?;
            } else {
                stage.file.file().write_all(bytes).map_err(io)?;
            }
            stage.input_hasher.update(bytes);
            stage.received += bytes.len() as u64;
            Ok(stage.received)
        })();
        match result {
            Ok(next) => {
                self.imports.insert(id.into(), stage);
                Ok(next)
            }
            Err(e) => Err(e),
        }
    }
    pub fn finish(&mut self, id: &str, max_blocks: usize) -> Result<NativeImportProgress, Error> {
        if max_blocks == 0 || max_blocks > 64 {
            return Err(bad("verification step requires 1..64 blocks"));
        }
        if let Some(dataset) = self.completed.get(id).and_then(|id| self.datasets.get(id)) {
            return Ok(done(id, &dataset.info));
        }
        let stage = self
            .imports
            .get(id)
            .ok_or_else(|| bad("unknown reference import"))?;
        if stage.received != stage.expected {
            return Err(bad(
                "input is not complete; append the remaining exact bytes",
            ));
        }
        if self.datasets.len() >= MAX_DATASETS {
            return Err(bad(
                "at most four reference datasets; delete one or cancel the import",
            ));
        }
        let mut stage = self.imports.remove(id).unwrap();
        let result = (|| {
            if stage.verification.is_none() {
                prepare(&mut stage)?;
            }
            let verification = stage.verification.as_mut().unwrap();
            for _ in 0..max_blocks {
                if verification.verified == verification.total {
                    break;
                }
                let contig = &verification.reference.manifest.contigs[verification.contig];
                let descriptor = &contig.blocks[verification.block];
                let bytes = read_at(
                    stage.file.file(),
                    descriptor.offset,
                    u64::from(descriptor.count).div_ceil(2) as usize,
                )?;
                let block = verification
                    .verifier
                    .verify(&contig.id, descriptor.index, &bytes)?;
                verification
                    .semantic
                    .as_mut()
                    .ok_or_else(|| bad("missing semantic verification state"))?
                    .push(&block)?;
                verification.verified += 1;
                verification.block += 1;
                if verification.block == contig.blocks.len() {
                    verification.contig += 1;
                    verification.block = 0;
                }
            }
            Ok(verification.verified == verification.total)
        })();
        match result {
            Err(e) => Err(e), // OwnedFile closes/deletes; no dataset can be partially installed.
            Ok(false) => {
                let verification = stage.verification.as_ref().unwrap();
                let progress = NativeImportProgress {
                    done: false,
                    import_id: id.into(),
                    source_bytes: stage.received,
                    blocks_verified: verification.verified,
                    blocks_total: verification.total,
                    dataset: None,
                };
                self.imports.insert(id.into(), stage);
                Ok(progress)
            }
            Ok(true) => {
                let mut verification = stage.verification.take().unwrap();
                let receipt = verification
                    .semantic
                    .take()
                    .ok_or_else(|| bad("missing semantic EOF validator"))?
                    .finish()?;
                if receipt.manifest_hash != verification.reference.manifest_hash
                    || receipt.content_revision != verification.reference.manifest.content_revision
                {
                    return Err(bad("semantic receipt identity mismatch"));
                }
                self.serial = self.serial.checked_add(1).ok_or(Error::TooLarge)?;
                let dataset_id = format!("reference-dataset-{}-{}", self.namespace, self.serial);
                let info = NativeDataset {
                    id: dataset_id.clone(),
                    reference: verification.reference,
                    semantic_verified: true,
                    input_bytes: stage.received,
                };
                let progress = done(id, &info);
                self.datasets.insert(
                    dataset_id.clone(),
                    Dataset {
                        info,
                        file: stage.file,
                    },
                );
                self.completed.insert(id.into(), dataset_id);
                Ok(progress)
            }
        }
    }
    pub fn inspect(
        &self,
        id: &str,
        manifest_hash: &str,
        content_revision: &str,
    ) -> Result<NativeDataset, Error> {
        Ok(self
            .dataset(id, manifest_hash, content_revision)?
            .info
            .clone())
    }
    fn dataset(
        &self,
        id: &str,
        manifest_hash: &str,
        content_revision: &str,
    ) -> Result<&Dataset, Error> {
        let d = self
            .datasets
            .get(id)
            .ok_or_else(|| bad("unknown reference dataset"))?;
        if d.info.reference.manifest_hash != manifest_hash
            || d.info.reference.manifest.content_revision != content_revision
        {
            return Err(bad("stale dataset manifest/content identity"));
        }
        Ok(d)
    }
    pub fn list(&self) -> Vec<NativeDataset> {
        self.datasets.values().map(|d| d.info.clone()).collect()
    }
    pub fn read_blocks(
        &mut self,
        id: &str,
        manifest_hash: &str,
        content_revision: &str,
        requests: &[BlockRequest],
    ) -> Result<Vec<EncodedBlock>, Error> {
        self.dataset(id, manifest_hash, content_revision)?;
        if requests.is_empty() || requests.len() > 64 {
            return Err(bad("request 1..64 encoded blocks"));
        }
        let dataset = self.datasets.get_mut(id).unwrap();
        let mut selected = vec![];
        for r in requests {
            let contig = dataset
                .info
                .reference
                .manifest
                .contigs
                .iter()
                .find(|c| c.id == r.contig_id)
                .ok_or_else(|| bad("unknown contig"))?;
            let descriptor = usize::try_from(r.index)
                .ok()
                .and_then(|i| contig.blocks.get(i))
                .ok_or_else(|| bad("unknown block"))?;
            selected.push((r.contig_id.clone(), descriptor.clone()));
        }
        let mut out = vec![];
        for (contig_id, descriptor) in selected {
            let encoded = read_at(
                dataset.file.file(),
                descriptor.offset,
                u64::from(descriptor.count).div_ceil(2) as usize,
            )?;
            if blake3::hash(&encoded).to_hex().as_str() != descriptor.encoded_hash {
                return Err(bad("managed block changed; current search cannot continue"));
            }
            out.push(EncodedBlock {
                contig_id,
                descriptor,
                encoded,
            });
        }
        Ok(out)
    }
    pub fn read_export(
        &mut self,
        id: &str,
        manifest_hash: &str,
        content_revision: &str,
        offset: u64,
        max_bytes: usize,
    ) -> Result<Vec<u8>, Error> {
        self.dataset(id, manifest_hash, content_revision)?;
        let dataset = self.datasets.get_mut(id).unwrap();
        if max_bytes == 0 || max_bytes > MAX_CHUNK || offset > dataset.info.reference.file_length {
            return Err(bad("export requires in-range offset and 1..1MiB bytes"));
        }
        let count = (dataset.info.reference.file_length - offset).min(max_bytes as u64) as usize;
        read_at(dataset.file.file(), offset, count)
    }
    pub fn cancel_import(&mut self, id: &str) -> Result<(), Error> {
        self.imports
            .remove(id)
            .ok_or_else(|| bad("unknown active reference import"))?;
        Ok(())
    }
    pub fn delete_dataset(
        &mut self,
        id: &str,
        manifest_hash: &str,
        content_revision: &str,
    ) -> Result<(), Error> {
        self.dataset(id, manifest_hash, content_revision)?;
        self.datasets.remove(id);
        self.completed.retain(|_, value| value != id);
        Ok(())
    }
}
impl Drop for NativeGenomes {
    fn drop(&mut self) {
        self.imports.clear();
        self.datasets.clear();
        // ManagedSession removes its empty child after all owned file handles close.
        // No absolute-path or recursive deletion is used here.
    }
}
fn done(id: &str, dataset: &NativeDataset) -> NativeImportProgress {
    let blocks = dataset
        .reference
        .manifest
        .contigs
        .iter()
        .map(|c| c.blocks.len() as u64)
        .sum();
    NativeImportProgress {
        done: true,
        import_id: id.into(),
        source_bytes: dataset.input_bytes,
        blocks_verified: blocks,
        blocks_total: blocks,
        dataset: Some(dataset.clone()),
    }
}
fn read_at(file: &mut File, offset: u64, count: usize) -> Result<Vec<u8>, Error> {
    file.seek(SeekFrom::Start(offset)).map_err(io)?;
    let mut bytes = vec![0; count];
    file.read_exact(&mut bytes).map_err(io)?;
    Ok(bytes)
}
fn write_events(file: &mut OwnedFile, events: &[EncodedBlock]) -> Result<(), Error> {
    for event in events {
        if file.file().metadata().map_err(io)?.len() != event.descriptor.offset {
            return Err(bad("append-only block offset mismatch"));
        }
        file.file().seek(SeekFrom::End(0)).map_err(io)?;
        file.file().write_all(&event.encoded).map_err(io)?;
    }
    Ok(())
}
fn prepare(stage: &mut Import) -> Result<(), Error> {
    if stage
        .expected_hash
        .as_ref()
        .is_some_and(|h| stage.input_hasher.finalize().to_hex().as_str() != h)
    {
        return Err(bad("uploaded source hash mismatch"));
    }
    let reference = if stage.format == ImportFormat::Fasta {
        let builder = stage
            .fasta
            .take()
            .ok_or_else(|| bad("missing streaming FASTA parser"))?;
        let (last, finished) = builder.finish()?;
        write_events(&mut stage.file, &last)?;
        stage
            .file
            .file()
            .write_all(&genome::encode_manifest(&finished.manifest)?)
            .map_err(io)?;
        stage
            .file
            .file()
            .write_all(&genome::encode_tail(&finished)?)
            .map_err(io)?;
        stage.file.file().flush().map_err(io)?;
        if stage.file.file().metadata().map_err(io)?.len() != finished.file_length {
            return Err(bad("generated reference file length mismatch"));
        }
        finished
    } else {
        stage.file.file().flush().map_err(io)?;
        let length = stage.file.file().metadata().map_err(io)?.len();
        if length != stage.expected
            || length < genome::HEADER.len() as u64 + genome::TAIL_BYTES as u64
        {
            return Err(bad("native reference file length mismatch"));
        }
        let header = read_at(stage.file.file(), 0, genome::HEADER.len())?;
        let tail = read_at(
            stage.file.file(),
            length - genome::TAIL_BYTES as u64,
            genome::TAIL_BYTES,
        )?;
        let tail = genome::parse_tail(&header, &tail, length)?;
        let bytes = read_at(
            stage.file.file(),
            tail.manifest_offset,
            usize::try_from(tail.manifest_length).map_err(|_| Error::TooLarge)?,
        )?;
        genome::decode_manifest(&bytes, &tail)?
    };
    if reference.file_length > stage.reserved {
        return Err(Error::TooLarge);
    }
    let total = reference
        .manifest
        .contigs
        .iter()
        .map(|c| c.blocks.len() as u64)
        .sum();
    stage.verification = Some(Verification {
        verifier: BlockVerifier::begin(&reference.manifest)?,
        semantic: Some(SemanticValidator::begin(&reference.manifest)?),
        reference,
        contig: 0,
        block: 0,
        verified: 0,
        total,
    });
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::{
        BlockRequest, BlockVerifier, BuildInfo, ImportFormat, MAX_CHUNK, NativeDataset,
        NativeGenomes, fs, genome,
    };
    fn finish(manager: &mut NativeGenomes, id: &str) -> NativeDataset {
        loop {
            let progress = manager.finish(id, 64).unwrap();
            if progress.done {
                return progress.dataset.unwrap();
            }
            assert!(progress.dataset.is_none());
            assert!(progress.blocks_verified < progress.blocks_total);
        }
    }
    fn small(manager: &mut NativeGenomes) -> NativeDataset {
        let bytes = b">a fixture\nACGTRYN\n";
        let id = manager
            .begin(
                ImportFormat::Fasta,
                bytes.len() as u64,
                Some(blake3::hash(bytes).to_hex().to_string()),
                None,
                None,
            )
            .unwrap();
        for (index, chunk) in bytes.chunks(3).enumerate() {
            manager.append(&id, (index * 3) as u64, chunk).unwrap();
        }
        finish(manager, &id)
    }
    #[test]
    fn offsets_capacity_and_poison_cleanup_do_not_install_partial_datasets() {
        let mut manager = NativeGenomes::new().unwrap();
        let a = manager
            .begin(ImportFormat::Fasta, 4, None, None, None)
            .unwrap();
        let b = manager
            .begin(ImportFormat::Fasta, 4, None, None, None)
            .unwrap();
        assert!(
            manager
                .begin(ImportFormat::Fasta, 4, None, None, None)
                .is_err()
        );
        assert!(manager.append(&a, 1, b">a\nA").is_err());
        assert_eq!(manager.append(&a, 0, b">a\nA").unwrap(), 4);
        manager.cancel_import(&b).unwrap();
        let dataset = finish(&mut manager, &a);
        assert!(dataset.semantic_verified);
        assert_eq!(
            manager.finish(&a, 1).unwrap().dataset.unwrap().id,
            dataset.id
        );
        assert!(
            manager
                .inspect(
                    &dataset.id,
                    &"0".repeat(64),
                    &dataset.reference.manifest.content_revision
                )
                .is_err()
        );
        let invalid = manager
            .begin(ImportFormat::Fasta, 4, None, None, None)
            .unwrap();
        assert!(manager.append(&invalid, 0, b"bad!").is_err());
        assert!(manager.finish(&invalid, 1).is_err());
        assert_eq!(manager.list().len(), 1);
    }
    #[test]
    fn handles_are_bound_to_their_manager_epoch_even_for_identical_dna() {
        let mut first = NativeGenomes::new().unwrap();
        let mut second = NativeGenomes::new().unwrap();
        let a = first
            .begin(ImportFormat::Fasta, 4, None, None, None)
            .unwrap();
        let b = second
            .begin(ImportFormat::Fasta, 4, None, None, None)
            .unwrap();
        assert_ne!(a, b);
        assert!(second.append(&a, 0, b">a\nA").is_err());
        second.append(&b, 0, b">a\nA").unwrap();
        first.cancel_import(&a).unwrap();
        second.cancel_import(&b).unwrap();
        let a = small(&mut first);
        let b = small(&mut second);
        assert_eq!(a.reference.manifest_hash, b.reference.manifest_hash);
        assert_ne!(a.id, b.id);
        assert!(
            second
                .inspect(
                    &a.id,
                    &a.reference.manifest_hash,
                    &a.reference.manifest.content_revision
                )
                .is_err()
        );
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            assert_eq!(
                fs::metadata(&first.directory).unwrap().permissions().mode() & 0o777,
                0o700
            );
        }
    }
    #[test]
    fn source_hash_and_native_block_corruption_fail_before_publication() {
        let mut manager = NativeGenomes::new().unwrap();
        let raw = b">a\nACGT\n";
        let bad = manager
            .begin(
                ImportFormat::Fasta,
                raw.len() as u64,
                Some("0".repeat(64)),
                None,
                None,
            )
            .unwrap();
        manager.append(&bad, 0, raw).unwrap();
        assert!(manager.finish(&bad, 1).is_err());
        assert!(manager.list().is_empty());
        let dataset = small(&mut manager);
        let mut bytes = manager
            .read_export(
                &dataset.id,
                &dataset.reference.manifest_hash,
                &dataset.reference.manifest.content_revision,
                0,
                MAX_CHUNK,
            )
            .unwrap();
        bytes[genome::HEADER.len()] ^= 1;
        let bad = manager
            .begin(ImportFormat::Cgref, bytes.len() as u64, None, None, None)
            .unwrap();
        manager.append(&bad, 0, &bytes).unwrap();
        assert!(manager.finish(&bad, 64).is_err());
        assert_eq!(manager.list().len(), 1);
        assert!(
            manager
                .begin(
                    ImportFormat::Cgref,
                    100,
                    None,
                    Some(BuildInfo::default()),
                    None
                )
                .is_err()
        );
    }
    #[test]
    fn deleting_and_dropping_close_files_before_removing_owned_storage() {
        let mut manager = NativeGenomes::new().unwrap();
        let directory = manager.directory.clone();
        let dataset = small(&mut manager);
        assert!(directory.exists());
        manager
            .delete_dataset(
                &dataset.id,
                &dataset.reference.manifest_hash,
                &dataset.reference.manifest.content_revision,
            )
            .unwrap();
        assert_eq!(fs::read_dir(&directory).unwrap().count(), 0);
        let pending = manager
            .begin(ImportFormat::Fasta, 100, None, None, None)
            .unwrap();
        manager.append(&pending, 0, b">a\nACGT").unwrap();
        drop(manager);
        assert!(!directory.exists());
    }
    #[test]
    fn streamed_seventy_mib_fasta_and_over_thirty_two_mib_native_file_never_need_whole_file_ram() {
        let mut manager = NativeGenomes::new().unwrap();
        let directory = manager.directory.clone();
        let header = b">large_contig streamed-unwrapped\n";
        let bases = 70u64 * 1024 * 1024 + 3;
        let id = manager
            .begin(
                ImportFormat::Fasta,
                bases + header.len() as u64,
                None,
                None,
                None,
            )
            .unwrap();
        manager.append(&id, 0, header).unwrap();
        let chunk = vec![b'A'; MAX_CHUNK];
        let mut supplied = 0u64;
        while supplied < bases {
            let n = (bases - supplied).min(chunk.len() as u64) as usize;
            manager
                .append(&id, header.len() as u64 + supplied, &chunk[..n])
                .unwrap();
            supplied += n as u64;
        }
        let first = finish(&mut manager, &id);
        assert_eq!(first.reference.manifest.total_bases, bases);
        assert!(first.reference.file_length > 32 * 1024 * 1024);
        assert!(first.reference.manifest.contigs[0].length > 10_000_000);
        let native = manager
            .begin(
                ImportFormat::Cgref,
                first.reference.file_length,
                None,
                None,
                None,
            )
            .unwrap();
        let mut offset = 0u64;
        while offset < first.reference.file_length {
            let bytes = manager
                .read_export(
                    &first.id,
                    &first.reference.manifest_hash,
                    &first.reference.manifest.content_revision,
                    offset,
                    MAX_CHUNK,
                )
                .unwrap();
            assert!(bytes.len() <= MAX_CHUNK);
            manager.append(&native, offset, &bytes).unwrap();
            offset += bytes.len() as u64;
        }
        let second = finish(&mut manager, &native);
        assert_eq!(first.reference, second.reference);
        assert!(second.semantic_verified);
        let contig = &second.reference.manifest.contigs[0];
        let last = contig.blocks.last().unwrap();
        let blocks = manager
            .read_blocks(
                &second.id,
                &second.reference.manifest_hash,
                &second.reference.manifest.content_revision,
                &[BlockRequest {
                    contig_id: contig.id.clone(),
                    index: last.index,
                }],
            )
            .unwrap();
        let verified = BlockVerifier::begin(&second.reference.manifest)
            .unwrap()
            .verify(&contig.id, last.index, &blocks[0].encoded)
            .unwrap();
        assert_eq!(verified.dna(), b"AAA");
        assert!(
            manager
                .read_export(
                    &second.id,
                    &second.reference.manifest_hash,
                    &second.reference.manifest.content_revision,
                    second.reference.file_length,
                    1
                )
                .unwrap()
                .is_empty()
        );
        drop(manager);
        assert!(!directory.exists());
    }
}
