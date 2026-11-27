//! Cheap editable references to bounded documentary recipes, never live material.
use collagene_core::{
    run_lineage::EventBinding,
    run_material::{Capsule, MaterialOrigin, MaterialSelection},
};
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::sync::Arc;

#[derive(Clone, Debug)]
pub struct PinnedOrigin {
    material: Arc<MaterialOrigin>,
    revision: String,
}
impl PinnedOrigin {
    pub(crate) fn from_material(origin: MaterialOrigin) -> Result<Self, String> {
        origin.validate_documentary().map_err(|e| e.to_string())?;
        struct Hash {
            inner: blake3::Hasher,
            bytes: usize,
        }
        impl std::io::Write for Hash {
            fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
                self.bytes = self
                    .bytes
                    .checked_add(bytes.len())
                    .filter(|n| *n <= 128 * 1024 * 1024)
                    .ok_or_else(|| {
                        std::io::Error::other("Origin identity exceeds bounded source context")
                    })?;
                self.inner.update(bytes);
                Ok(bytes.len())
            }
            fn flush(&mut self) -> std::io::Result<()> {
                Ok(())
            }
        }
        let mut hash = Hash {
            inner: blake3::Hasher::new(),
            bytes: 0,
        };
        serde_json::to_writer(&mut hash, &origin).map_err(|e| e.to_string())?;
        Ok(Self {
            material: Arc::new(origin),
            revision: hash.inner.finalize().to_hex().to_string(),
        })
    }
    /// Immutable source identity computed once when captured or decoded; no replay.
    pub(crate) fn revision(&self) -> &str {
        &self.revision
    }
    pub(crate) fn as_ref(&self) -> &MaterialOrigin {
        &self.material
    }
    pub(crate) fn event(&self) -> &EventBinding {
        self.material.event()
    }
    pub(crate) fn selection(&self) -> MaterialSelection {
        self.material.selection()
    }
    pub(crate) fn validate_documentary(&self) -> Result<(), String> {
        self.material
            .validate_documentary()
            .map_err(|e| e.to_string())
    }
    /// Advanced run-only workflows must refuse derived digest material explicitly.
    pub(crate) fn as_run(&self) -> Result<&Capsule, String> {
        match self.as_ref() {
            MaterialOrigin::Run(origin) => Ok(origin),
            MaterialOrigin::Derived(_) => Err("This workflow does not yet support target-free derived products. Use ordinary PCR or an assembly preparation to retain their actual material origin.".into()),
        }
    }
    pub(crate) fn label(&self) -> String {
        match self.as_ref() {
            MaterialOrigin::Run(p) => format!(
                "Run {}… · product {}",
                p.selection.run_id.chars().take(12).collect::<String>(),
                p.selection.step_index + 1
            ),
            MaterialOrigin::Derived(_) => format!(
                "Derived material · event {}…",
                self.event().event_id.chars().take(12).collect::<String>()
            ),
        }
    }
}
impl PartialEq for PinnedOrigin {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.material, &other.material) || self.revision == other.revision
    }
}
impl Eq for PinnedOrigin {}
impl Serialize for PinnedOrigin {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.as_ref().serialize(serializer)
    }
}
impl<'de> Deserialize<'de> for PinnedOrigin {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        let origin = MaterialOrigin::deserialize(deserializer)?;
        Self::from_material(origin).map_err(serde::de::Error::custom)
    }
}
