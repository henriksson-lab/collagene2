//! Small, versioned display preferences, separate from workspace history or DNA.
use serde::{Deserialize, Serialize};

pub const MAX_BYTES: usize = 1024;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ViewMode {
    Split,
    Sequence,
    Map,
}

impl ViewMode {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Split => "split",
            Self::Sequence => "sequence",
            Self::Map => "map",
        }
    }

    /// Optional parsing for an explicit display control, never a DOM authority.
    #[allow(clippy::should_implement_trait)]
    pub fn from_str(value: &str) -> Option<Self> {
        match value {
            "split" => Some(Self::Split),
            "sequence" => Some(Self::Sequence),
            "map" => Some(Self::Map),
            _ => None,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Preferences {
    pub version: u8,
    pub library_width: u16,
    pub split_per_mille: u16,
    pub view: ViewMode,
    pub library_visible: bool,
}

impl Default for Preferences {
    fn default() -> Self {
        Self {
            version: 1,
            library_width: 244,
            split_per_mille: 524,
            view: ViewMode::Split,
            library_visible: true,
        }
    }
}

fn validate(preferences: &Preferences) -> Result<(), String> {
    if preferences.version != 1 {
        return Err("Unsupported workspace view preference version".into());
    }
    if !(180..=360).contains(&preferences.library_width) {
        return Err("Workspace library width must be between 180 and 360".into());
    }
    if !(100..=900).contains(&preferences.split_per_mille) {
        return Err("Workspace pane ratio must be between 100 and 900 per mille".into());
    }
    Ok(())
}

pub fn decode(bytes: &[u8]) -> Result<Preferences, String> {
    if bytes.len() > MAX_BYTES {
        return Err("Workspace view preferences exceed 1024 bytes".into());
    }
    let preferences = serde_json::from_slice(bytes)
        .map_err(|error| format!("Invalid workspace view preferences: {error}"))?;
    validate(&preferences)?;
    Ok(preferences)
}

pub fn encode(preferences: &Preferences) -> Result<Vec<u8>, String> {
    validate(preferences)?;
    let bytes = serde_json::to_vec(preferences)
        .map_err(|error| format!("Cannot encode workspace view preferences: {error}"))?;
    if bytes.len() > MAX_BYTES {
        return Err("Workspace view preferences exceed 1024 bytes".into());
    }
    Ok(bytes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn literal_defaults_and_all_modes_round_trip_at_supported_boundaries() {
        let literal = br#"{"version":1,"library_width":244,"split_per_mille":524,"view":"split","library_visible":true}"#;
        assert_eq!(decode(literal).unwrap(), Preferences::default());
        assert_eq!(encode(&Preferences::default()).unwrap(), literal);
        for (view, width, ratio, visible) in [
            (ViewMode::Split, 180, 100, false),
            (ViewMode::Sequence, 360, 900, true),
            (ViewMode::Map, 244, 524, false),
        ] {
            let preferences = Preferences {
                version: 1,
                library_width: width,
                split_per_mille: ratio,
                view,
                library_visible: visible,
            };
            assert_eq!(ViewMode::from_str(view.as_str()), Some(view));
            let bytes = encode(&preferences).unwrap();
            assert!(bytes.len() <= MAX_BYTES);
            assert_eq!(decode(&bytes).unwrap(), preferences);
        }
        assert_eq!(ViewMode::from_str("Map"), None);
        assert_eq!(ViewMode::from_str(" map "), None);
    }

    #[test]
    fn unsupported_versions_and_outside_bounds_are_rejected_on_read_and_write() {
        for preferences in [
            Preferences {
                version: 0,
                ..Preferences::default()
            },
            Preferences {
                version: 2,
                ..Preferences::default()
            },
            Preferences {
                library_width: 179,
                ..Preferences::default()
            },
            Preferences {
                library_width: 361,
                ..Preferences::default()
            },
            Preferences {
                split_per_mille: 99,
                ..Preferences::default()
            },
            Preferences {
                split_per_mille: 901,
                ..Preferences::default()
            },
        ] {
            assert!(encode(&preferences).is_err());
            // Independently serialized malformed state must not bypass read validation.
            assert!(decode(&serde_json::to_vec(&preferences).unwrap()).is_err());
        }
    }

    #[test]
    fn unknown_fields_types_missing_fields_and_duplicate_fields_do_not_restore() {
        for literal in [
            r#"{"version":1,"library_width":244,"split_per_mille":524,"view":"split","library_visible":true,"extra":0}"#,
            r#"{"version":1,"library_width":244,"split_per_mille":524,"view":"Split","library_visible":true}"#,
            r#"{"version":1,"library_width":-1,"split_per_mille":524,"view":"map","library_visible":true}"#,
            r#"{"version":1,"library_width":244,"split_per_mille":524,"view":"map","library_visible":"true"}"#,
            r#"{"version":1,"library_width":244,"split_per_mille":524,"view":"map"}"#,
            r#"{"version":1,"version":1,"library_width":244,"split_per_mille":524,"view":"map","library_visible":true}"#,
        ] {
            assert!(decode(literal.as_bytes()).is_err(), "accepted {literal}");
        }
    }

    #[test]
    fn oversized_invalid_bytes_fail_before_json_decoding() {
        assert_eq!(
            decode(&vec![0xff; MAX_BYTES + 1]).unwrap_err(),
            "Workspace view preferences exceed 1024 bytes"
        );
        assert!(
            decode(&[0xff])
                .unwrap_err()
                .contains("Invalid workspace view preferences")
        );
    }
}
