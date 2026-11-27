//! Stable recipe identities; document order and display labels are not identity.
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Role {
    pub input_id: String,
    pub role: String,
}
fn token(value: &str) -> bool {
    !value.is_empty()
        && value.len() <= 128
        && value
            .bytes()
            .all(|v| v.is_ascii_alphanumeric() || matches!(v, b'-' | b'_' | b'.' | b':'))
}
pub fn resolve(value: &Option<Vec<Role>>, count: usize) -> Result<Vec<Role>, String> {
    if count == 0 || count > 12 {
        return Err("Choose 1–12 ordered assembly sources.".into());
    }
    let Some(roles) = value else {
        return Ok((0..count)
            .map(|index| Role {
                input_id: format!("source-{}", index + 1),
                role: format!("fragment-{}", index + 1),
            })
            .collect());
    };
    if roles.len() != count {
        return Err("Stable assembly identities do not match the ordered source slots.".into());
    }
    let mut inputs = std::collections::BTreeSet::new();
    let mut fragments = std::collections::BTreeSet::new();
    for role in roles {
        if !token(&role.input_id)
            || !token(&role.role)
            || !inputs.insert(&role.input_id)
            || !fragments.insert(&role.role)
        {
            return Err("Assembly source and fragment identities must be valid and unique.".into());
        }
    }
    Ok(roles.clone())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn ordinary_positional_roles_are_explicit_and_named_roles_follow_their_sources() {
        let old = resolve(&None, 2).unwrap();
        assert_eq!(
            old[0],
            Role {
                input_id: "source-1".into(),
                role: "fragment-1".into()
            }
        );
        let mut stable = vec![
            Role {
                input_id: "choice-promoter-a".into(),
                role: "bin-promoter".into(),
            },
            Role {
                input_id: "choice-cds-z".into(),
                role: "bin-cds".into(),
            },
        ];
        stable.swap(0, 1);
        let result = resolve(&Some(stable), 2).unwrap();
        assert_eq!(result[0].input_id, "choice-cds-z");
        assert_eq!(result[1].role, "bin-promoter");
    }
    #[test]
    fn mismatched_duplicate_and_invalid_identities_are_not_reassigned_by_index() {
        let role = Role {
            input_id: "source-a".into(),
            role: "bin-a".into(),
        };
        assert!(resolve(&Some(vec![role.clone()]), 2).is_err());
        assert!(resolve(&Some(vec![role.clone(), role.clone()]), 2).is_err());
        let mut other = role.clone();
        other.role = "bin-b".into();
        assert!(resolve(&Some(vec![role.clone(), other]), 2).is_err());
        let mut invalid = role;
        invalid.input_id = "source with spaces".into();
        assert!(resolve(&Some(vec![invalid]), 1).is_err());
        assert!(resolve(&None, 0).is_err());
        assert!(resolve(&None, 13).is_err());
    }
}
