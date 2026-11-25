//! Portable editable pool settings, never a saved inspection or verification result.
use collagene_core::golden_gate::{GoldenGateEnzyme, GoldenGateOptions};
use serde::{Deserialize, Deserializer, Serialize, de};
use std::fmt;

pub const MAX_CHOICES: usize = 12;
pub const MAX_FIELD_BYTES: usize = 20;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PoolForm {
    pub reaction: bool,
    /// Per-slot substrate: 0 prepared, 1 whole original, 2 whole amplicon.
    /// Empty is the legacy/default form; initialize from the actual slot count.
    #[serde(deserialize_with = "deserialize_choices")]
    pub choices: Vec<u8>,
    pub circular: bool,
    pub both: bool,
    #[serde(deserialize_with = "deserialize_field")]
    pub minimum: String,
    #[serde(deserialize_with = "deserialize_field")]
    pub maximum: String,
    #[serde(deserialize_with = "deserialize_field")]
    pub products: String,
    #[serde(deserialize_with = "deserialize_field")]
    pub states: String,
}

impl Default for PoolForm {
    fn default() -> Self {
        Self {
            reaction: false,
            choices: Vec::new(),
            circular: true,
            both: true,
            minimum: "1".into(),
            maximum: "16".into(),
            products: "128".into(),
            states: "50000".into(),
        }
    }
}

fn valid_field(value: &str) -> bool {
    value.len() <= MAX_FIELD_BYTES && !value.chars().any(char::is_control)
}

impl PoolForm {
    /// Validate safe editable state. Empty or temporarily invalid numeric text
    /// remains valid draft data; it cannot be compiled into a runnable request.
    pub fn validate(&self) -> Result<(), String> {
        if self.choices.len() > MAX_CHOICES || self.choices.iter().any(|&choice| choice > 2) {
            return Err("Pool settings require at most 12 known substrate choices.".into());
        }
        if [&self.minimum, &self.maximum, &self.products, &self.states]
            .into_iter()
            .any(|value| !valid_field(value))
        {
            return Err(
                "Pool limit fields allow at most 20 bytes and no control characters.".into(),
            );
        }
        Ok(())
    }

    /// Parse a runnable options draft. These are the existing native finite-pool
    /// bounds, not proof that substrates, PCR choices, targets or chemistry are
    /// valid. The caller must still use prepared/reaction native preflight on
    /// the complete source-bound request before inspection.
    pub fn options(&self, enzyme: GoldenGateEnzyme) -> Result<GoldenGateOptions, String> {
        self.validate()?;
        let options = GoldenGateOptions {
            enzyme,
            circular: self.circular,
            both_orientations: self.both,
            min_fragments: self
                .minimum
                .parse()
                .map_err(|_| "Enter a minimum fragment count.")?,
            max_fragments: self
                .maximum
                .parse()
                .map_err(|_| "Enter a maximum fragment count.")?,
            max_products: self
                .products
                .parse()
                .map_err(|_| "Enter a product limit.")?,
            max_states: self
                .states
                .parse()
                .map_err(|_| "Enter a search-state limit.")?,
        };
        // GoldenGateOptions has no public validation method. Keep this UI
        // parse guard aligned with golden_gate_reaction::preflight and
        // prepared_golden_gate::preflight; those remain authoritative.
        if options.min_fragments == 0
            || options.min_fragments > options.max_fragments
            || options.max_fragments > 32
            || options.max_products == 0
            || options.max_products > 1000
            || options.max_states == 0
            || options.max_states > 1_000_000
        {
            return Err("Pool limits require fragments 1–32, products 1–1000 and search states 1–1,000,000.".into());
        }
        Ok(options)
    }
}

fn deserialize_choices<'de, D: Deserializer<'de>>(deserializer: D) -> Result<Vec<u8>, D::Error> {
    struct ChoicesVisitor;
    impl<'de> de::Visitor<'de> for ChoicesVisitor {
        type Value = Vec<u8>;
        fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
            f.write_str("at most 12 substrate choices from 0, 1, 2")
        }
        fn visit_seq<A: de::SeqAccess<'de>>(self, mut seq: A) -> Result<Self::Value, A::Error> {
            if seq.size_hint().is_some_and(|len| len > MAX_CHOICES) {
                return Err(de::Error::custom("too many pool substrate choices"));
            }
            let mut choices = Vec::with_capacity(2);
            while choices.len() < MAX_CHOICES {
                let Some(choice) = seq.next_element::<u8>()? else {
                    return Ok(choices);
                };
                if choice > 2 {
                    return Err(de::Error::custom("unknown pool substrate choice"));
                }
                choices.push(choice);
            }
            // JSON sequences have no size hint. Decode at most one extra byte,
            // never collect an unbounded sequence or arbitrary ignored value.
            if seq.next_element::<u8>()?.is_some() {
                return Err(de::Error::custom("too many pool substrate choices"));
            }
            Ok(choices)
        }
    }
    deserializer.deserialize_seq(ChoicesVisitor)
}

fn deserialize_field<'de, D: Deserializer<'de>>(deserializer: D) -> Result<String, D::Error> {
    struct FieldVisitor;
    impl de::Visitor<'_> for FieldVisitor {
        type Value = String;
        fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
            f.write_str("a pool limit draft of at most 20 bytes without control characters")
        }
        fn visit_str<E: de::Error>(self, value: &str) -> Result<String, E> {
            if !valid_field(value) {
                return Err(E::custom(
                    "pool limit field exceeds bounds or contains controls",
                ));
            }
            Ok(value.into())
        }
        fn visit_string<E: de::Error>(self, value: String) -> Result<String, E> {
            if !valid_field(&value) {
                return Err(E::custom(
                    "pool limit field exceeds bounds or contains controls",
                ));
            }
            Ok(value)
        }
    }
    // Avoid requesting an owned String before checking its size. Formats may
    // still allocate their own scratch space (e.g. escaped JSON); enclosing
    // archive/message byte limits must bound the serialized input as usual.
    deserializer.deserialize_str(FieldVisitor)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn decode(value: serde_json::Value) -> Result<PoolForm, serde_json::Error> {
        serde_json::from_str(&value.to_string())
    }
    fn value() -> serde_json::Value {
        serde_json::to_value(PoolForm::default()).unwrap()
    }

    #[test]
    fn default_matches_native_options_without_retaining_evidence() {
        let form = PoolForm::default();
        assert!(form.choices.is_empty());
        assert!(!form.reaction);
        assert_eq!(
            form.options(GoldenGateEnzyme::BsaI).unwrap(),
            GoldenGateOptions::default()
        );
        let restored: PoolForm =
            serde_json::from_str(&serde_json::to_string(&form).unwrap()).unwrap();
        assert_eq!(restored, form);
    }

    #[test]
    fn editable_incomplete_text_roundtrips_but_cannot_run() {
        for incomplete in ["", "-", "not ready", "1.", " 1"] {
            let mut form = PoolForm::default();
            form.minimum = incomplete.into();
            form.validate().unwrap();
            assert_eq!(decode(serde_json::to_value(&form).unwrap()).unwrap(), form);
            assert!(form.options(GoldenGateEnzyme::BsaI).is_err());
        }
    }

    #[test]
    fn restored_settings_keep_mode_geometry_and_slot_material_choices() {
        let mut form = PoolForm::default();
        form.reaction = true;
        form.circular = false;
        form.both = false;
        form.choices = vec![2, 0, 1];
        let restored = decode(serde_json::to_value(&form).unwrap()).unwrap();
        assert_eq!(restored, form);
        let options = restored.options(GoldenGateEnzyme::BsmBI).unwrap();
        assert_eq!(options.enzyme, GoldenGateEnzyme::BsmBI);
        assert!(!options.circular && !options.both_orientations);
    }

    #[test]
    fn unknown_evidence_fields_are_rejected() {
        for key in [
            "report",
            "verified",
            "fresh",
            "candidate_index",
            "assessment",
        ] {
            let mut json = value();
            json[key] = serde_json::json!(true);
            assert!(decode(json).is_err());
        }
    }

    #[test]
    fn json_choice_sequences_are_bounded_without_a_size_hint() {
        let mut json = value();
        json["choices"] = serde_json::json!(vec![2; MAX_CHOICES]);
        assert_eq!(decode(json.clone()).unwrap().choices.len(), MAX_CHOICES);
        json["choices"] = serde_json::json!(vec![0; MAX_CHOICES + 1]);
        assert!(decode(json).is_err());
        for choices in [
            serde_json::json!([3]),
            serde_json::json!([-1]),
            serde_json::json!([256]),
            serde_json::json!(["1"]),
        ] {
            let mut json = value();
            json["choices"] = choices;
            assert!(decode(json).is_err());
        }
    }

    #[test]
    fn numeric_text_has_byte_and_control_bounds_on_each_field() {
        for key in ["minimum", "maximum", "products", "states"] {
            let mut json = value();
            json[key] = serde_json::json!("0".repeat(MAX_FIELD_BYTES));
            assert!(decode(json.clone()).is_ok());
            for text in [
                "0".repeat(MAX_FIELD_BYTES + 1),
                "é".repeat(11),
                "1\n".into(),
                "\u{7f}".into(),
                "\u{0085}".into(),
            ] {
                json[key] = serde_json::json!(text);
                assert!(decode(json.clone()).is_err());
            }
        }
        let mut form = PoolForm::default();
        form.states = "\t1".into();
        assert!(form.validate().is_err());
        assert!(form.options(GoldenGateEnzyme::BsaI).is_err());
    }

    #[test]
    fn numeric_limits_do_not_turn_an_editable_draft_into_runnable_options() {
        for (field, text) in [
            ("minimum", "0"),
            ("minimum", "17"),
            ("maximum", "0"),
            ("maximum", "33"),
            ("products", "0"),
            ("products", "1001"),
            ("states", "0"),
            ("states", "1000001"),
            ("states", "99999999999999999999"),
        ] {
            let mut json = value();
            json[field] = serde_json::json!(text);
            let form = decode(json).unwrap();
            form.validate().unwrap();
            assert!(
                form.options(GoldenGateEnzyme::BsaI).is_err(),
                "{field}={text}"
            );
        }
        let mut form = PoolForm::default();
        form.minimum = "32".into();
        form.maximum = "32".into();
        form.products = "1000".into();
        form.states = "1000000".into();
        assert!(form.options(GoldenGateEnzyme::BsaI).is_ok());
    }

    #[test]
    fn runnable_ui_limits_agree_with_native_reaction_preflight() {
        use collagene_core::{Molecule, Topology, golden_gate_reaction as reaction};
        let source = Molecule::dna("source", "ACGT", Topology::Linear).unwrap();
        let target = Molecule::dna("target", "ACGT", Topology::Linear).unwrap();
        for (field, number) in [
            ("minimum", 0),
            ("minimum", 1),
            ("minimum", 17),
            ("maximum", 0),
            ("maximum", 32),
            ("maximum", 33),
            ("products", 0),
            ("products", 1000),
            ("products", 1001),
            ("states", 0),
            ("states", 1_000_000),
            ("states", 1_000_001),
        ] {
            let mut json = value();
            json[field] = serde_json::json!(number.to_string());
            let form = decode(json).unwrap();
            // Construct the corresponding raw scientific request independently
            // of PoolForm::options so this catches future native-bound drift.
            let mut native = GoldenGateOptions::default();
            match field {
                "minimum" => native.min_fragments = number,
                "maximum" => native.max_fragments = number,
                "products" => native.max_products = number,
                "states" => native.max_states = number,
                _ => unreachable!(),
            }
            let request = reaction::Request {
                name: "native limits".into(),
                rows: vec![reaction::Row {
                    role: "substrate".into(),
                    input_id: "source".into(),
                    substrate: reaction::Substrate::WholeSource,
                }],
                options: native,
            };
            assert_eq!(
                form.options(GoldenGateEnzyme::BsaI).is_ok(),
                reaction::preflight(&[("source", &source)], &target, &request).is_ok(),
                "{field}={number}",
            );
        }
    }
}
