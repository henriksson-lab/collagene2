//! Documentary download of current recorded physical/reference evidence. Reading
//! attachments never reconstructs material or establishes current verification.
use crate::{run_lineage, workspace::Workspace};
use collagene_model::Error;
use serde::{Deserialize, Serialize};
use std::io::Write;
pub const MAX_BYTES: usize = 32 * 1024 * 1024;
#[derive(Clone, Debug, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Download {
    pub run_id: String,
    pub recorded_engine: String,
    pub step_count: usize,
    pub bytes: Option<Vec<u8>>,
}
struct Bounded(Vec<u8>);
impl Write for Bounded {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if bytes.len() > MAX_BYTES.saturating_sub(4096).saturating_sub(self.0.len()) {
            return Err(std::io::Error::other(
                "scientific run evidence exceeds download limit; export the native workspace",
            ));
        }
        self.0
            .try_reserve(bytes.len())
            .map_err(|_| std::io::Error::other("scientific evidence allocation"))?;
        self.0.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
pub fn download(workspace: &Workspace, run_id: &str) -> Result<Download, Error> {
    let record = run_lineage::load_run(workspace, run_id)?;
    let physical = crate::physical_run_evidence::find(workspace, run_id)?;
    let reference = crate::reference_run_evidence::find(workspace, run_id)?;
    let step_count = physical.as_ref().map_or(0, |r| r.steps.len())
        + reference.as_ref().map_or(0, |r| r.steps.len());
    let bytes = if step_count == 0 {
        None
    } else {
        let mut output = Bounded(vec![]);
        output
            .write_all(b"{\"format\":\"scientific-run-bundle-v1\",\"run_id\":")
            .map_err(|_| Error::TooLarge)?;
        serde_json::to_writer(&mut output, run_id).map_err(|_| Error::TooLarge)?;
        output
            .write_all(b",\"attachments\":[")
            .map_err(|_| Error::TooLarge)?;
        let mut count = 0;
        for (hash, attachment) in &workspace.attachments {
            if attachment.kind != crate::physical_run_evidence::KIND
                && attachment.kind != crate::reference_run_evidence::KIND
            {
                continue;
            }
            #[derive(Deserialize)]
            struct Owner {
                run_id: String,
            }
            let owner: Owner = serde_json::from_slice(&attachment.bytes)
                .map_err(|e| Error::Format(e.to_string()))?;
            if owner.run_id != run_id {
                continue;
            }
            if count > 0 {
                output.write_all(b",").map_err(|_| Error::TooLarge)?;
            }
            output
                .write_all(b"{\"kind\":")
                .map_err(|_| Error::TooLarge)?;
            serde_json::to_writer(&mut output, &attachment.kind).map_err(|_| Error::TooLarge)?;
            output
                .write_all(b",\"blake3\":")
                .map_err(|_| Error::TooLarge)?;
            serde_json::to_writer(&mut output, hash).map_err(|_| Error::TooLarge)?;
            output
                .write_all(b",\"evidence\":")
                .map_err(|_| Error::TooLarge)?;
            output
                .write_all(&attachment.bytes)
                .map_err(|_| Error::TooLarge)?;
            output.write_all(b"}").map_err(|_| Error::TooLarge)?;
            count += 1;
        }
        output.write_all(b"]}").map_err(|_| Error::TooLarge)?;
        Some(output.0)
    };
    Ok(Download {
        run_id: run_id.into(),
        recorded_engine: record.assessment.engine_version,
        step_count,
        bytes,
    })
}
