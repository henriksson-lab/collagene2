//! Cheap editable references to bounded documentary recipes, never live material.
use collagene_core::run_material::Capsule;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::{ops::Deref, sync::Arc};

#[derive(Clone, Debug)]
pub struct PinnedOrigin(Arc<Capsule>);
impl PinnedOrigin {
    pub(crate) fn new(origin: Capsule) -> Result<Self, String> {
        origin.validate_documentary().map_err(|e| e.to_string())?;
        Ok(Self(Arc::new(origin)))
    }
    pub(crate) fn as_ref(&self) -> &Capsule {
        &self.0
    }
}
impl Deref for PinnedOrigin {
    type Target = Capsule;
    fn deref(&self) -> &Capsule {
        &self.0
    }
}
impl PartialEq for PinnedOrigin {
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.0, &other.0) || self.0 == other.0
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
        Self::new(Capsule::deserialize(deserializer)?).map_err(serde::de::Error::custom)
    }
}
