//! Operator-owned individual credentials. No role or author claim comes from wire.
use crate::{AccessKey, Error, Result};
use collagene_lab_library::collaboration::{ActorId, Author, Role};
use serde::Deserialize;
use std::{
    path::{Path, PathBuf},
    sync::{Arc, RwLock},
};
pub const MAX_MEMBERS: usize = 256;
const MAX_CONFIG_BYTES: u64 = 256 * 1024;
pub struct MemberCredential {
    key: AccessKey,
    author: Author,
    role: Role,
}
impl MemberCredential {
    pub fn new(key: AccessKey, author: Author, role: Role) -> Result<Self> {
        let author = Author::new(ActorId::parse(author.id.as_str())?, author.display_name)?;
        Ok(Self { key, author, role })
    }
}
#[derive(Debug, Clone)]
pub struct AuthenticatedMember {
    pub author: Author,
    pub role: Role,
}
impl AuthenticatedMember {
    pub fn permits_library(&self, request: &collagene_library_protocol::Request) -> bool {
        self.role != Role::Reader || !request.is_mutation()
    }
    pub fn permits_sync(&self, request: &collagene_library_sync::Request) -> bool {
        self.role != Role::Reader
            || matches!(request, collagene_library_sync::Request::Lookup { .. })
    }
}
/// Clones share the current registry. A trusted host can atomically replace it;
/// accepted requests resolve current membership again immediately before dispatch.
#[derive(Clone)]
pub struct MemberRegistry(Arc<RwLock<Vec<MemberCredential>>>);
impl MemberRegistry {
    pub fn new(members: Vec<MemberCredential>) -> Result<Self> {
        Self::validate(&members)?;
        Ok(Self(Arc::new(RwLock::new(members))))
    }
    fn validate(members: &[MemberCredential]) -> Result<()> {
        if members.is_empty() || members.len() > MAX_MEMBERS {
            return Err(Error::Configuration("member count"));
        }
        for (i, m) in members.iter().enumerate() {
            for previous in &members[..i] {
                if m.author.id == previous.author.id || m.key.same_secret(&previous.key) {
                    return Err(Error::Configuration(
                        "duplicate member identity or credential",
                    ));
                }
            }
        }
        Ok(())
    }
    pub fn replace(&self, members: Vec<MemberCredential>) -> Result<()> {
        Self::validate(&members)?;
        *self
            .0
            .write()
            .map_err(|_| Error::Configuration("member registry unavailable"))? = members;
        Ok(())
    }
    pub fn authenticate(&self, authorization: &[u8]) -> Option<AuthenticatedMember> {
        let members = self.0.read().ok()?;
        let mut found = None;
        // Evaluate every key; do not reveal registry position by early return.
        for m in members.iter() {
            if m.key.authenticates(authorization) {
                found = Some(AuthenticatedMember {
                    author: m.author.clone(),
                    role: m.role,
                });
            }
        }
        found
    }
    /// Identity-only revocation check for bounded hosted-session cleanup. This
    /// never substitutes for credential authentication on an incoming request.
    pub(crate) fn contains_actor(&self, actor: &str) -> bool {
        self.0.read().is_ok_and(|members| {
            members
                .iter()
                .any(|member| member.author.id.as_str() == actor)
        })
    }
    /// Reads a bounded owner-only registry and individually protected key files.
    /// Relative key paths are resolved against the registry's own directory.
    pub fn read_file(path: &Path, assets_root: &Path) -> Result<Self> {
        use cap_fs_ext::{FollowSymlinks, OpenOptionsFollowExt};
        use cap_std::{
            ambient_authority,
            fs::{Dir, OpenOptions},
        };
        use std::io::Read;
        let assets_root = assets_root.canonicalize()?;
        if path.canonicalize()?.starts_with(&assets_root) {
            return Err(Error::Configuration(
                "member registry must be outside static assets",
            ));
        }
        let parent = path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let directory = Dir::open_ambient_dir(parent, ambient_authority())?;
        let mut options = OpenOptions::new();
        options.read(true).follow(FollowSymlinks::No);
        let mut file = directory
            .open_with(
                path.file_name()
                    .ok_or(Error::Configuration("member registry path"))?,
                &options,
            )?
            .into_std();
        let metadata = file.metadata()?;
        if !metadata.is_file() || metadata.len() > MAX_CONFIG_BYTES {
            return Err(Error::Configuration("member registry size/type"));
        }
        #[cfg(unix)]
        {
            use std::os::unix::fs::MetadataExt;
            if metadata.mode() & 0o077 != 0 {
                return Err(Error::Configuration("member registry must be owner-only"));
            }
        }
        let mut bytes = Vec::new();
        Read::by_ref(&mut file)
            .take(MAX_CONFIG_BYTES + 1)
            .read_to_end(&mut bytes)?;
        if bytes.len() as u64 > MAX_CONFIG_BYTES {
            return Err(Error::Configuration("member registry size"));
        }
        let config: Config = serde_json::from_slice(&bytes)
            .map_err(|_| Error::Configuration("member registry fields"))?;
        if config.schema_version != 1
            || config.members.is_empty()
            || config.members.len() > MAX_MEMBERS
        {
            return Err(Error::Configuration("member registry version/count"));
        }
        let mut members = Vec::new();
        for m in config.members {
            let key_path = if m.key_file.is_absolute() {
                m.key_file
            } else {
                parent.join(m.key_file)
            };
            if key_path.canonicalize()?.starts_with(&assets_root) {
                return Err(Error::Configuration(
                    "member keys must be outside static assets",
                ));
            }
            members.push(MemberCredential::new(
                AccessKey::read(&key_path)?,
                Author::new(ActorId::parse(m.actor_id)?, m.display_name)?,
                m.role,
            )?);
        }
        Self::new(members)
    }
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Config {
    schema_version: u16,
    members: Vec<ConfigMember>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct ConfigMember {
    actor_id: String,
    display_name: String,
    role: Role,
    key_file: PathBuf,
}
