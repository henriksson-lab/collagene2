//! Explicit native remote retirement. Local biological data is never deleted.
use crate::{Error, RemoteConfig, Repository, Result, durable, state_write, store, wire};
use serde::{Deserialize, Serialize};
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RemovalDescription {
    pub local_lab: wire::LabId,
    pub remote: RemoteConfig,
    pub expected_checkpoint: Option<store::Digest>,
    pub document_bindings: u32,
    pub folder_bindings: u32,
    pub completed_tree_jobs: u16,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RemovalReceipt {
    pub schema_version: u16,
    pub local_lab: wire::LabId,
    pub remote: RemoteConfig,
    pub review_digest: wire::Digest,
    pub removed_document_bindings: u32,
    pub removed_folder_bindings: u32,
    pub removed_completed_tree_jobs: u16,
}
#[derive(Debug, Clone)]
struct JobPin {
    key: durable::CheckpointKey,
    digest: store::Digest,
}
#[derive(Debug, Clone)]
pub struct RemovalReview {
    description: RemovalDescription,
    digest: wire::Digest,
    jobs: Vec<JobPin>,
    index_digest: Option<store::Digest>,
    remaining_jobs: Vec<String>,
}
#[derive(Debug, Clone)]
pub struct ApprovedRemoteRemoval(RemovalReview);
impl RemovalReview {
    pub fn description(&self) -> &RemovalDescription {
        &self.description
    }
    pub fn digest(&self) -> wire::Digest {
        self.digest.clone()
    }
    pub fn approve(&self) -> ApprovedRemoteRemoval {
        ApprovedRemoteRemoval(self.clone())
    }
}
fn receipt_key(remote: &str) -> Result<durable::CheckpointKey> {
    if remote.len() != 36
        || !remote.starts_with("rem-")
        || !remote[4..]
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    {
        return Err(Error::Invalid("remote ID"));
    }
    Ok(durable::CheckpointKey::repository(format!(
        "remote-removal-{remote}"
    )))
}
fn job_key(id: &str) -> Result<durable::CheckpointKey> {
    if id.len() != 37
        || !id.starts_with("tree-")
        || !id[5..]
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    {
        return Err(Error::Invalid("tree job ID"));
    }
    Ok(durable::CheckpointKey::repository(format!("job-{id}")))
}
impl Repository {
    /// Secret-free durable outcome. It conveys local bookkeeping completion only,
    /// never successful vault deletion, remote data deletion or scientific status.
    pub fn remote_removal_receipt(&mut self, remote: &str) -> Result<Option<RemovalReceipt>> {
        let key = receipt_key(remote)?;
        let Some(c) = self.db.read_checkpoint(&self.access, &key)? else {
            return Ok(None);
        };
        if c.value.version != 1 || c.value.bytes.len() > 16 * 1024 {
            return Err(Error::Invalid("remote removal receipt"));
        }
        let r: RemovalReceipt = serde_json::from_slice(&c.value.bytes)
            .map_err(|_| Error::Invalid("remote removal receipt"))?;
        r.remote.validate()?;
        if r.schema_version != 1
            || r.remote.id != remote
            || r.local_lab != self.local_context()?.lab_id
            || r.removed_document_bindings > 100_000
            || r.removed_folder_bindings > 100_000
            || r.removed_completed_tree_jobs > 128
        {
            return Err(Error::Invalid("remote removal receipt scope"));
        }
        Ok(Some(r))
    }
    /// Capture exact native scope. No remote network/vault/working-file operations.
    pub fn prepare_remote_removal(&mut self, remote: &str) -> Result<RemovalReview> {
        self.reload()?;
        let config = self
            .state
            .remotes
            .get(remote)
            .cloned()
            .ok_or(Error::Invalid("unknown remote"))?;
        if self
            .state
            .pending
            .values()
            .any(|pending| pending.remote.id == remote)
        {
            return Err(Error::Pending);
        }
        if self.remote_removal_receipt(remote)?.is_some() {
            return Err(Error::Invalid("removed remote identity was reused"));
        }
        let docs: Vec<_> = self
            .state
            .documents
            .iter()
            .filter_map(|(id, bindings)| bindings.get(remote).map(|binding| (id, binding)))
            .collect();
        let folders: Vec<_> = self
            .state
            .folders
            .iter()
            .filter_map(|(id, bindings)| bindings.get(remote).map(|binding| (id, binding)))
            .collect();
        let binding_digest = wire::Digest::of(
            &serde_json::to_vec(&(&docs, &folders))
                .map_err(|_| Error::Invalid("remote binding snapshot"))?,
        );
        let document_bindings = docs
            .len()
            .try_into()
            .map_err(|_| Error::Invalid("binding count"))?;
        let folder_bindings = folders
            .len()
            .try_into()
            .map_err(|_| Error::Invalid("binding count"))?;
        let index = self.db.read_checkpoint(
            &self.access,
            &durable::CheckpointKey::repository("tree-jobs"),
        )?;
        let ids: Vec<String> = match &index {
            None => vec![],
            Some(c) => {
                if c.value.version != 1 || c.value.bytes.len() > 16384 {
                    return Err(Error::Invalid("tree job index"));
                }
                serde_json::from_slice(&c.value.bytes)
                    .map_err(|_| Error::Invalid("tree job index"))?
            }
        };
        if ids.len() > 128 {
            return Err(Error::Invalid("tree job count"));
        }
        let mut seen = std::collections::BTreeSet::new();
        let mut jobs = vec![];
        let mut remaining = vec![];
        for id in ids {
            let key = job_key(&id)?;
            if !seen.insert(id.clone()) {
                return Err(Error::Invalid("duplicate tree job"));
            }
            let job = self
                .tree_clone_job(&id)
                .map_err(|_| Error::Invalid("tree job snapshot"))?;
            let c = self
                .db
                .read_checkpoint(&self.access, &key)?
                .ok_or(Error::Invalid("missing tree job"))?;
            // Public loader enforces the complete bounded job contract. Verify its exact
            // bytes still match this later CAS pin rather than mixing concurrent reads.
            if c.value.version != 1
                || serde_json::to_vec(&job).map_err(|_| Error::Invalid("tree job snapshot"))?
                    != c.value.bytes
            {
                return Err(Error::Stale);
            }
            if job.remote.id == remote {
                if job.remote != config {
                    return Err(Error::Stale);
                }
                if job.phase != crate::folders::TreePhase::Complete {
                    return Err(Error::Invalid(
                        "finish the selected remote's tree clones before removal",
                    ));
                }
                jobs.push(JobPin {
                    key,
                    digest: c.digest,
                });
            } else {
                remaining.push(id);
            }
        }
        let description = RemovalDescription {
            local_lab: self.local_context()?.lab_id,
            remote: config,
            expected_checkpoint: self.checkpoint.clone(),
            document_bindings,
            folder_bindings,
            completed_tree_jobs: jobs
                .len()
                .try_into()
                .map_err(|_| Error::Invalid("tree job count"))?,
        };
        let pins: Vec<_> = jobs.iter().map(|j| (&j.key.key, &j.digest)).collect();
        let index_digest = index.map(|c| c.digest);
        let digest = wire::Digest::of(
            &serde_json::to_vec(&(
                "collagene-remote-removal-1",
                &description,
                binding_digest,
                &index_digest,
                &pins,
                &remaining,
            ))
            .map_err(|_| Error::Invalid("remote removal binding"))?,
        );
        Ok(RemovalReview {
            description,
            digest,
            jobs,
            index_digest,
            remaining_jobs: remaining,
        })
    }
    /// Host orders any credential cleanup separately using the captured full remote
    /// scope. Never claim a vault and SQLite share one atomic transaction.
    pub fn remove_remote_reviewed(
        &mut self,
        approved: &ApprovedRemoteRemoval,
    ) -> Result<RemovalReceipt> {
        let review = &approved.0;
        let remote = &review.description.remote.id;
        if let Some(receipt) = self.remote_removal_receipt(remote)? {
            if receipt.review_digest == review.digest
                && receipt.local_lab == review.description.local_lab
            {
                self.reload()?;
                return Ok(receipt);
            }
            return Err(Error::Stale);
        }
        let current = self.prepare_remote_removal(remote)?;
        if current.digest != review.digest {
            return Err(Error::Stale);
        }
        let mut next = self.state.clone();
        next.remotes.remove(remote);
        for bindings in next.documents.values_mut() {
            bindings.remove(remote);
        }
        next.documents.retain(|_, bindings| !bindings.is_empty());
        for bindings in next.folders.values_mut() {
            bindings.remove(remote);
        }
        next.folders.retain(|_, bindings| !bindings.is_empty());
        let receipt = RemovalReceipt {
            schema_version: 1,
            local_lab: review.description.local_lab.clone(),
            remote: review.description.remote.clone(),
            review_digest: review.digest.clone(),
            removed_document_bindings: review.description.document_bindings,
            removed_folder_bindings: review.description.folder_bindings,
            removed_completed_tree_jobs: review.description.completed_tree_jobs,
        };
        let mut writes = vec![state_write(
            &next,
            current.description.expected_checkpoint.clone(),
        )?];
        if !current.jobs.is_empty() {
            writes.push(durable::CheckpointWrite {
                key: durable::CheckpointKey::repository("tree-jobs"),
                expected: current.index_digest.clone(),
                value: Some(durable::CheckpointValue {
                    version: 1,
                    bytes: serde_json::to_vec(&current.remaining_jobs)
                        .map_err(|_| Error::Invalid("tree job index"))?,
                }),
            });
            for job in current.jobs {
                writes.push(durable::CheckpointWrite {
                    key: job.key,
                    expected: Some(job.digest),
                    value: None,
                });
            }
        } else if current.index_digest.is_some() {
            // Assert the index still matches even when it contributes no deletion. A new
            // selected-remote job cannot race removal without its own source CAS failure.
            writes.push(durable::CheckpointWrite {
                key: durable::CheckpointKey::repository("tree-jobs"),
                expected: current.index_digest.clone(),
                value: Some(durable::CheckpointValue {
                    version: 1,
                    bytes: serde_json::to_vec(&current.remaining_jobs)
                        .map_err(|_| Error::Invalid("tree job index"))?,
                }),
            });
        } else {
            // Materialize the empty index, asserting absence as well. Tree-job
            // creators must also assert the repository state checkpoint in the
            // same transaction so a retired remote cannot acquire a late job.
            writes.push(durable::CheckpointWrite {
                key: durable::CheckpointKey::repository("tree-jobs"),
                expected: None,
                value: Some(durable::CheckpointValue {
                    version: 1,
                    bytes: b"[]".to_vec(),
                }),
            });
        }
        writes.push(durable::CheckpointWrite {
            key: receipt_key(remote)?,
            expected: None,
            value: Some(durable::CheckpointValue {
                version: 1,
                bytes: serde_json::to_vec(&receipt)
                    .map_err(|_| Error::Invalid("removal receipt"))?,
            }),
        });
        self.db
            .commit_checkpoint_retirements(&self.access, &writes)?;
        self.reload()?;
        Ok(receipt)
    }
}
