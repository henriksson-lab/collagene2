//! Persisted declarations only; profile acceptance and component selection are not proof.
use collagene_core::{
    multisite_recombination as science,
    recombination::{ProfileSnapshot, reference},
};
use serde::{Deserialize, Serialize};
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ComponentIntent {
    pub batch_revision: String,
    pub candidate_id: usize,
    pub component_id: usize,
    pub source_revision: String,
}
#[derive(Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Form {
    pub profile_json: String,
    pub reaction: reference::Reaction,
    pub specificities: String,
    pub reviewed: bool,
    pub selection: Option<ComponentIntent>,
}
impl Default for Form {
    fn default() -> Self {
        let profile_json = ProfileSnapshot::manufacturer_context().ok().and_then(|s| {
            let p=s.profile;
            serde_json::to_string_pretty(&reference::multisite::Profile {id:p.id,edition:p.edition,sources:p.sources,variants:p.variants,limitations:vec!["Registered deposited contexts for specificities1 and2 only; not a complete MultiSite kit or a physical activity profile. Import a reviewed full-context profile for other specificities.".into()]}).ok()
        }).unwrap_or_default();
        Self {
            profile_json,
            reaction: reference::Reaction::LR,
            specificities: "1,2".into(),
            reviewed: false,
            selection: None,
        }
    }
}
fn digest(s: &str) -> bool {
    s.len() == 64 && s.bytes().all(|b| b.is_ascii_hexdigit())
}
impl Form {
    pub fn three_entry() -> Result<Self, String> {
        let profile =
            reference::multisite::three_entry_source_contexts().map_err(|e| e.to_string())?;
        Ok(Self {
            profile_json: serde_json::to_string_pretty(&profile).map_err(|e| e.to_string())?,
            reaction: reference::Reaction::LR,
            specificities: "4,1,2,3".into(),
            reviewed: false,
            selection: None,
        })
    }
    pub fn validate_shape(&self) -> Result<(), String> {
        if self.profile_json.len() > 512 * 1024 || self.specificities.len() > 128 {
            return Err("Reference profile exceeds its bounded editable fields.".into());
        }
        if let Some(s) = &self.selection {
            if !digest(&s.batch_revision)
                || !digest(&s.source_revision)
                || s.candidate_id >= 256
                || s.component_id >= 32
            {
                return Err("Invalid reference component intent.".into());
            }
        }
        Ok(())
    }
    pub fn request(&self) -> Result<science::Request, String> {
        self.validate_shape()?;
        if !self.reviewed {
            return Err("Review and acknowledge this exact full-context profile first.".into());
        }
        let profile: reference::multisite::Profile = serde_json::from_str(&self.profile_json)
            .map_err(|e| format!("Invalid imported reference profile: {e}"))?;
        let specificities = self
            .specificities
            .split(',')
            .map(|s| {
                s.trim()
                    .parse::<u8>()
                    .map_err(|_| "Use comma-separated specificity IDs1..16.".to_string())
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut limits = reference::Limits::default();
        limits.max_components_per_candidate = 32;
        let request = science::Request {
            reaction: reference::multisite::Request {
                profile,
                reaction: self.reaction,
                specificities,
                limits,
            },
            annotation_limits: science::AnnotationLimits::default(),
        };
        science::request_preflight(&request).map_err(|e| e.to_string())?;
        Ok(request)
    }
    pub fn selection(&self) -> Option<collagene_core::combinatorial_design::MultiSiteSelection> {
        self.selection.as_ref().map(
            |s| collagene_core::combinatorial_design::MultiSiteSelection {
                selection: science::Selection {
                    batch_revision: s.batch_revision.clone(),
                    candidate_id: s.candidate_id,
                    component_id: s.component_id,
                },
                source_revision: s.source_revision.clone(),
            },
        )
    }
}
pub const ORDINARY_ID: &str = "ordinary-multi";
pub fn ordinary_roles(count: usize) -> Result<Vec<crate::assembly_role_ids::Role>, String> {
    if !(1..=12).contains(&count) {
        return Err("Editable assembly identities require1..12 source slots.".into());
    }
    (0..count)
        .map(|i| {
            let (input_id, role) = collagene_core::combinatorial_design::source_role_ids(
                ORDINARY_ID,
                &format!("{ORDINARY_ID}-row-{}", count * 2 + 1),
                &format!("{ORDINARY_ID}-bin-{}", i * 2 + 1),
                &format!("{ORDINARY_ID}-choice-{}", i * 2 + 2),
            )
            .map_err(|e| e.to_string())?;
            Ok(crate::assembly_role_ids::Role { input_id, role })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn presets_require_explicit_review_and_pin_the_declared_contexts() {
        let mut standard = Form::default();
        assert!(standard.request().is_err());
        standard.reviewed = true;
        let request = standard.request().unwrap();
        assert_eq!(request.reaction.specificities, vec![1, 2]);
        assert_eq!(request.reaction.limits.max_components_per_candidate, 32);
        let mut three = Form::three_entry().unwrap();
        assert!(!three.reviewed);
        assert!(three.selection.is_none());
        assert!(three.request().is_err());
        three.reviewed = true;
        let request = three.request().unwrap();
        assert_eq!(request.reaction.reaction, reference::Reaction::LR);
        assert_eq!(request.reaction.specificities, vec![4, 1, 2, 3]);
        assert_eq!(request.reaction.profile.variants.len(), 9);
        assert!(
            request
                .reaction
                .profile
                .variants
                .iter()
                .any(|v| v.family == reference::Family::P
                    && v.subtype == 3
                    && v.sequence.len() == 232)
        );
        assert!(
            request
                .reaction
                .profile
                .sources
                .iter()
                .any(|s| s.vector_name == "p5E-MCS")
        );
    }
}
