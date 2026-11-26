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

pub(crate) fn apply(
    state: &crate::cloning_methods_model::MethodsState,
    index: usize,
    document: usize,
    event: &collagene_core::run_lineage::EventBinding,
    inspection: &collagene_core::run_material::Inspection,
    revision: &str,
) -> Result<crate::cloning_methods_model::MethodsState, String> {
    use crate::cloning_methods_model::OverlapPreparation;
    let row = state
        .rows
        .get(index)
        .ok_or("The preparation row is unavailable.")?;
    if row.document != document
        || matches!(
            row.preparation,
            OverlapPreparation::DeclaredDuplex | OverlapPreparation::AnnealedOligos
        )
    {
        return Err("The source or preparation changed; select the saved material again.".into());
    }
    inspection
        .validate_binding(revision)
        .map_err(|e| e.to_string())?;
    if &inspection.event != event {
        return Err("The reconstruction belongs to a different source event.".into());
    }
    let origin = inspection
        .selected
        .as_ref()
        .ok_or("This originating run has no freshly reconstructable material.")?;
    if &origin.selection.event != event {
        return Err("The originating material differs from the exact selected event.".into());
    }
    let pin = PinnedOrigin::new(origin.clone())?;
    let mut next = state.clone();
    next.rows[index].run_origin = Some(pin);
    Ok(next)
}
