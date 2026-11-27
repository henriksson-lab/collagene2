//! Cheap editable references to bounded documentary recipes, never live material.
use collagene_core::{
    run_lineage::EventBinding,
    run_material::{Capsule, MaterialOrigin, MaterialSelection},
};
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::sync::Arc;

#[derive(Clone, Debug)]
pub struct PinnedOrigin(Arc<MaterialOrigin>);
impl PinnedOrigin {
    pub(crate) fn from_material(origin: MaterialOrigin) -> Result<Self, String> {
        origin.validate_documentary().map_err(|e| e.to_string())?;
        Ok(Self(Arc::new(origin)))
    }
    pub(crate) fn as_ref(&self) -> &MaterialOrigin {
        &self.0
    }
    pub(crate) fn event(&self) -> &EventBinding {
        self.0.event()
    }
    pub(crate) fn selection(&self) -> MaterialSelection {
        self.0.selection()
    }
    pub(crate) fn validate_documentary(&self) -> Result<(), String> {
        self.0.validate_documentary().map_err(|e| e.to_string())
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
        Arc::ptr_eq(&self.0, &other.0)
            || serde_json::to_vec(self.as_ref())
                .ok()
                .zip(serde_json::to_vec(other.as_ref()).ok())
                .is_some_and(|(a, b)| a == b)
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
