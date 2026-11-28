//! Explicit original-source format declarations; no retained scientific authority.
use super::model::Draft;
use crate::cloning::preparation::cloning_methods_model::{MethodsState, OverlapPreparation};
use collagene_core::{
    cloning::assembly::{
        intent::Selector,
        standards::declaration::{Declaration, DeclaredPart, ProfileRef},
    },
    material::lineage::EventBinding,
};
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Row {
    pub role: String,
    pub input_id: String,
    pub event: EventBinding,
    pub selector: Selector,
    pub preparation: OverlapPreparation,
    pub preparation_reverse: bool,
    pub slot: String,
    pub body_start: String,
    pub body_end: String,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StandardDraft {
    pub profile: ProfileRef,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Row, 12>")]
    pub rows: Vec<Row>,
}
impl StandardDraft {
    pub fn validate_storage(&self) -> Result<(), String> {
        if self.rows.len() > 12
            || self.profile.profile_id.is_empty()
            || self.profile.profile_id.len() > 128
            || self.profile.profile_identity.len() != 64
        {
            return Err("Invalid bounded standard declaration.".into());
        }
        for row in &self.rows {
            if row.role.is_empty()
                || row.role.len() > 128
                || row.input_id.is_empty()
                || row.input_id.len() > 128
                || row.slot.len() > 96
                || row.body_start.len() > 20
                || row.body_end.len() > 20
                || row.event.document_id.len() > 128
                || row.event.event_id.len() > 80
                || row.event.revision.len() != 64
            {
                return Err(
                    "Standard source/body declarations exceed their storage bounds.".into(),
                );
            }
        }
        Ok(())
    }
    fn validate_binding(
        row: &Row,
        slot: &super::model::Slot,
        preparation: &crate::cloning::preparation::cloning_methods_model::OverlapInput,
    ) -> Result<(), String> {
        if slot.captured.as_ref().map(|c| &c.source.binding) != Some(&row.event)
            || slot.selector != row.selector
            || preparation.preparation != row.preparation
            || preparation.reverse != row.preparation_reverse
            || slot.material_origin.is_some()
            || slot.digest_geometry.is_some()
            || preparation.run_origin.is_some()
        {
            return Err("A standard source event, interval, orientation or preparation changed; explicitly rebind or clear its declaration.".into());
        }
        Ok(())
    }
    pub fn validate_choices(
        &self,
        library: &super::assembly_alternatives::Library,
        construct: &super::assembly_alternatives::Construct,
    ) -> Result<(), String> {
        self.validate_storage()?;
        let roles = construct.edits.methods.role_ids(library.bins.len())?;
        if roles.len() != self.rows.len() {
            return Err("Standard roles changed; rebind or clear the declaration.".into());
        }
        for (((bin, choice_id), role), preparation) in library
            .bins
            .iter()
            .zip(&construct.choices)
            .zip(roles)
            .zip(&construct.edits.methods.rows)
        {
            let choice = bin
                .choices
                .iter()
                .find(|choice| &choice.id == choice_id)
                .ok_or("Declared source choice is unavailable.")?;
            let row = self
                .rows
                .iter()
                .find(|row| row.role == role.role && row.input_id == role.input_id)
                .ok_or("Standard role identity changed; rebind or clear it.")?;
            Self::validate_binding(row, &choice.slot, preparation)?;
        }
        Ok(())
    }

    pub fn capture(draft: &Draft, profile: ProfileRef) -> Result<Self, String> {
        let roles = draft.methods.role_ids(draft.slots.len())?;
        let rows = draft.slots.iter().zip(&draft.methods.rows).zip(roles).map(|((slot, preparation), role)| {
            let captured = slot.captured.as_ref().ok_or("Capture every declared original source first.")?;
            slot.selector.reference_interval().map_err(|e| e.to_string())?;
            if slot.material_origin.is_some() || slot.digest_geometry.is_some() || preparation.run_origin.is_some() {
                return Err("Standard source-format review needs original reference coordinates; saved/material strand coordinates are unsupported.".into());
            }
            Ok(Row {
                role: role.role, input_id: role.input_id, event: captured.source.binding.clone(), selector: slot.selector.clone(),
                preparation: preparation.preparation, preparation_reverse: preparation.reverse,
                slot: String::new(), body_start: String::new(), body_end: String::new(),
            })
        }).collect::<Result<Vec<_>, String>>()?;
        Ok(Self { profile, rows })
    }
    /// Reordering retains exact identities, but source/geometry changes need an
    /// explicit new binding; old numeric body offsets are never reinterpreted.
    pub fn validate_sources(&self, draft: &Draft) -> Result<(), String> {
        let roles = draft.methods.role_ids(draft.slots.len())?;
        if self.rows.len() != roles.len() {
            return Err(
                "Declared standard roles changed; rebind the declaration or clear it.".into(),
            );
        }
        for ((slot, preparation), role) in draft.slots.iter().zip(&draft.methods.rows).zip(&roles) {
            let row = self
                .rows
                .iter()
                .find(|row| row.role == role.role && row.input_id == role.input_id)
                .ok_or("A declared standard role was removed/reassigned; rebind or clear it.")?;
            Self::validate_binding(row, slot, preparation)?;
        }
        Ok(())
    }
    pub fn declaration(&self, methods: &MethodsState) -> Result<Declaration, String> {
        let roles = methods.role_ids(methods.rows.len())?;
        if roles.len() != self.rows.len() || self.rows.len() > 12 {
            return Err("Standard declaration differs from current ordered roles.".into());
        }
        let parts = roles
            .iter()
            .zip(&methods.rows)
            .map(|(role, preparation)| {
                let row = self
                    .rows
                    .iter()
                    .find(|row| row.role == role.role && row.input_id == role.input_id)
                    .ok_or("A standard role binding is unavailable; rebind or clear it.")?;
                if row.preparation != preparation.preparation
                    || row.preparation_reverse != preparation.reverse
                    || preparation.run_origin.is_some()
                {
                    return Err(
                        "A preparation changed after the standard was bound; rebind or clear it."
                            .into(),
                    );
                }
                let (start, end, reverse) = row
                    .selector
                    .reference_interval()
                    .map_err(|e| e.to_string())?;
                fn number(value: &str) -> Result<usize, String> {
                    if value.is_empty()
                        || value.len() > 20
                        || !value.bytes().all(|b| b.is_ascii_digit())
                    {
                        return Err(
                            "Enter explicit zero-based body boundaries for every standard part."
                                .into(),
                        );
                    }
                    value
                        .parse()
                        .map_err(|_| "Body boundary exceeds the supported range.".into())
                }
                Ok(DeclaredPart {
                    role: row.role.clone(),
                    input_id: row.input_id.clone(),
                    source_revision: row.event.revision.clone(),
                    start,
                    length: end
                        .checked_sub(start)
                        .ok_or("Invalid declared source interval.")?,
                    reverse,
                    slot: row.slot.clone(),
                    body_start: number(&row.body_start)?,
                    body_end: number(&row.body_end)?,
                })
            })
            .collect::<Result<Vec<_>, String>>()?;
        let declaration = Declaration {
            profile: self.profile.clone(),
            parts,
            limits: Default::default(),
        };
        declaration.validate().map_err(|e| e.to_string())?;
        Ok(declaration)
    }
}

/// Exact documentary pins also cover inactive construct declarations, including
/// roles whose current source choices were later removed. No equal-DNA relabeling.
pub(crate) fn declared_events(draft: &Draft) -> Vec<&EventBinding> {
    std::iter::once(&draft.methods)
        .chain(
            draft
                .alternatives
                .iter()
                .flat_map(|library| library.rows.iter().map(|row| &row.edits.methods)),
        )
        .flat_map(|methods| {
            methods
                .standard
                .iter()
                .flat_map(|standard| standard.rows.iter().map(|row| &row.event))
        })
        .collect()
}
