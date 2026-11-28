//! Typed admission precedes any JSON tree; the second pass never builds an input tree.
use super::{Error, bad, format_error};
use serde::{
    Deserialize, Serialize,
    de::{self, DeserializeOwned, DeserializeSeed, MapAccess, SeqAccess, Visitor},
};
use serde_json::Value;
use std::{
    collections::{BTreeMap, BTreeSet},
    fmt,
};

pub(super) fn json<T: DeserializeOwned>(bytes: &[u8]) -> Result<T, Error> {
    let parsed = serde_json::from_slice(bytes).map_err(format_error)?;
    let mut reader = serde_json::Deserializer::from_slice(bytes);
    Unique.deserialize(&mut reader).map_err(format_error)?;
    reader.end().map_err(format_error)?;
    Ok(parsed)
}
pub(super) fn json_exact<T: DeserializeOwned + Serialize>(bytes: &[u8]) -> Result<T, Error> {
    let parsed: T = serde_json::from_slice(bytes).map_err(format_error)?;
    let expected = serde_json::to_value(&parsed).map_err(format_error)?;
    let mut reader = serde_json::Deserializer::from_slice(bytes);
    Equal(&expected)
        .deserialize(&mut reader)
        .map_err(|_| bad("workspace typed data contains unsupported or noncanonical fields"))?;
    reader.end().map_err(format_error)?;
    Ok(parsed)
}
struct Unique;
impl<'de> DeserializeSeed<'de> for Unique {
    type Value = ();
    fn deserialize<D: serde::Deserializer<'de>>(self, d: D) -> Result<(), D::Error> {
        d.deserialize_any(self)
    }
}
impl<'de> Visitor<'de> for Unique {
    type Value = ();
    fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str("JSON without duplicate keys")
    }
    fn visit_unit<E: de::Error>(self) -> Result<(), E> {
        Ok(())
    }
    fn visit_bool<E: de::Error>(self, _: bool) -> Result<(), E> {
        Ok(())
    }
    fn visit_i64<E: de::Error>(self, _: i64) -> Result<(), E> {
        Ok(())
    }
    fn visit_u64<E: de::Error>(self, _: u64) -> Result<(), E> {
        Ok(())
    }
    fn visit_f64<E: de::Error>(self, _: f64) -> Result<(), E> {
        Ok(())
    }
    fn visit_str<E: de::Error>(self, _: &str) -> Result<(), E> {
        Ok(())
    }
    fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<(), A::Error> {
        while seq.next_element_seed(Unique)?.is_some() {}
        Ok(())
    }
    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<(), A::Error> {
        let mut seen = BTreeSet::new();
        while let Some(key) = map.next_key::<String>()? {
            if !seen.insert(key) {
                return Err(de::Error::custom("duplicate JSON object key"));
            }
            map.next_value_seed(Unique)?;
        }
        Ok(())
    }
}
struct Equal<'a>(&'a Value);
impl<'de> DeserializeSeed<'de> for Equal<'_> {
    type Value = ();
    fn deserialize<D: serde::Deserializer<'de>>(self, d: D) -> Result<(), D::Error> {
        d.deserialize_any(self)
    }
}
impl Equal<'_> {
    fn scalar<E: de::Error>(self, value: Value) -> Result<(), E> {
        if self.0 == &value {
            Ok(())
        } else {
            Err(E::custom("noncanonical JSON value"))
        }
    }
}
struct End;
impl<'de> DeserializeSeed<'de> for End {
    type Value = ();
    fn deserialize<D: serde::Deserializer<'de>>(self, _: D) -> Result<(), D::Error> {
        Err(de::Error::custom("unexpected JSON array element"))
    }
}
impl<'de> Visitor<'de> for Equal<'_> {
    type Value = ();
    fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
        f.write_str("the exact typed JSON value")
    }
    fn visit_unit<E: de::Error>(self) -> Result<(), E> {
        self.scalar(Value::Null)
    }
    fn visit_bool<E: de::Error>(self, v: bool) -> Result<(), E> {
        self.scalar(v.into())
    }
    fn visit_i64<E: de::Error>(self, v: i64) -> Result<(), E> {
        self.scalar(v.into())
    }
    fn visit_u64<E: de::Error>(self, v: u64) -> Result<(), E> {
        self.scalar(v.into())
    }
    fn visit_f64<E: de::Error>(self, v: f64) -> Result<(), E> {
        self.scalar(
            serde_json::Number::from_f64(v)
                .ok_or_else(|| E::custom("non-finite JSON"))?
                .into(),
        )
    }
    fn visit_str<E: de::Error>(self, v: &str) -> Result<(), E> {
        if self.0.as_str() == Some(v) {
            Ok(())
        } else {
            Err(E::custom("noncanonical JSON text"))
        }
    }
    fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<(), A::Error> {
        let values = self
            .0
            .as_array()
            .ok_or_else(|| <A::Error as de::Error>::custom("unexpected JSON array"))?;
        for value in values {
            seq.next_element_seed(Equal(value))?
                .ok_or_else(|| <A::Error as de::Error>::custom("missing JSON array element"))?;
        }
        if seq.next_element_seed(End)?.is_some() {
            return Err(de::Error::custom("unexpected JSON array element"));
        }
        Ok(())
    }
    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<(), A::Error> {
        let values = self
            .0
            .as_object()
            .ok_or_else(|| <A::Error as de::Error>::custom("unexpected JSON object"))?;
        let mut seen = BTreeSet::new();
        while let Some(key) = map.next_key::<String>()? {
            let value = values
                .get(&key)
                .ok_or_else(|| <A::Error as de::Error>::custom("unexpected JSON field"))?;
            if !seen.insert(key) {
                return Err(de::Error::custom("duplicate JSON object key"));
            }
            map.next_value_seed(Equal(value))?;
        }
        if seen.len() != values.len() {
            return Err(de::Error::custom("missing JSON field"));
        }
        Ok(())
    }
}
/// Bound object cardinality before deserializing its next value. The wire map is unchanged.
pub(super) fn bounded_map<'de, D, T, const MAX: usize>(
    d: D,
) -> Result<BTreeMap<String, T>, D::Error>
where
    D: serde::Deserializer<'de>,
    T: Deserialize<'de>,
{
    struct Map<T, const MAX: usize>(std::marker::PhantomData<T>);
    impl<'de, T: Deserialize<'de>, const MAX: usize> Visitor<'de> for Map<T, MAX> {
        type Value = BTreeMap<String, T>;
        fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
            f.write_str("bounded archive map")
        }
        fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Self::Value, A::Error> {
            if map.size_hint().is_some_and(|n| n > MAX) {
                return Err(de::Error::custom("archive map resource budget"));
            }
            let mut values = BTreeMap::new();
            while let Some(key) = map.next_key::<String>()? {
                if values.len() == MAX {
                    return Err(de::Error::custom("archive map resource budget"));
                }
                if values.contains_key(&key) {
                    return Err(de::Error::custom("duplicate JSON object key"));
                }
                values.insert(key, map.next_value()?);
            }
            Ok(values)
        }
    }
    d.deserialize_map(Map::<T, MAX>(std::marker::PhantomData))
}
#[derive(Deserialize)]
pub(super) struct Header {
    pub format: Option<String>,
    pub version: Option<u64>,
    pub molecule_schema: Option<u64>,
    pub plan_schema: Option<u64>,
}

#[cfg(test)]
mod tests {
    use super::super::{FeatureList, StepList};
    use super::*;
    #[test]
    fn typed_step_bound_precedes_canonical_tree_and_extra_element_decode() {
        let operation = serde_json::to_string(&crate::planning::plan::Operation::Overlap {
            fragments: vec![],
            overlaps: vec![],
            circular: false,
        })
        .unwrap();
        let mut input = format!("[{}", vec![operation.as_str(); 1000].join(","));
        input.push_str(",{ ");
        // Inspect typed admission before the public archive boundary redacts
        // decoder details; EOF would mean the excess element was requested.
        let error = serde_json::from_slice::<StepList>(input.as_bytes())
            .err()
            .unwrap();
        assert!(error.is_data(), "{error}");
        assert!(error.to_string().contains("resource budget"), "{error}");
        assert!(json_exact::<StepList>(input.as_bytes()).is_err());
    }
    #[test]
    fn typed_feature_bound_precedes_extra_feature_construction() {
        let feature = serde_json::to_string(&collagene_model::Feature {
            kind: "x".into(),
            location: collagene_model::Location::range(0, 1),
            qualifiers: vec![],
            needs_review: false,
        })
        .unwrap();
        let mut input = format!("[{}", vec![feature.as_str(); 100_000].join(","));
        input.push_str(",{ ");
        let error = collagene_model::decode::with_native_frame_budget(
            super::super::MAX_EXPANDED_BYTES,
            || {
                serde_json::from_slice::<FeatureList>(input.as_bytes())
                    .err()
                    .unwrap()
            },
        );
        assert!(error.is_data(), "{error}");
        assert!(error.to_string().contains("resource budget"), "{error}");
        assert!(json_exact::<FeatureList>(input.as_bytes()).is_err());
    }
    #[test]
    fn archive_map_capacity_stops_before_the_next_value() {
        #[derive(Deserialize)]
        struct Record {
            #[serde(deserialize_with = "bounded_map::<_, String, 2>")]
            values: BTreeMap<String, String>,
        }
        let error = serde_json::from_slice::<Record>(br#"{"values":{"a":"x","b":"y","c":{ "#)
            .err()
            .unwrap();
        assert!(error.is_data(), "{error}");
        assert!(error.to_string().contains("resource budget"), "{error}");
        let record = serde_json::from_slice::<Record>(br#"{"values":{"a":"x","b":"y"}}"#).unwrap();
        assert_eq!(record.values.len(), 2);
    }
    #[test]
    fn metadata_record_content_is_typed_even_before_its_encoding_tag() {
        use super::super::ValueRecord;
        let error = serde_json::from_slice::<ValueRecord>(br#"{"blob":[ "#)
            .err()
            .unwrap();
        assert!(error.is_data(), "{error}");
        assert!(
            serde_json::from_slice::<ValueRecord>(br#"{"blob":"evidence","encoding":"utf8"}"#)
                .is_ok()
        );
        assert!(
            serde_json::from_slice::<ValueRecord>(br#"{"encoding":"base64","blob":"evidence"}"#)
                .is_ok()
        );
    }
    #[test]
    fn canonical_comparison_rejects_nested_duplicate_and_ignored_fields() {
        // Derive current location syntax rather than maintaining another encoding.
        let feature = collagene_model::Feature {
            kind: "x".into(),
            location: collagene_model::Location::range(0, 1),
            qualifiers: vec![],
            needs_review: false,
        };
        let current = serde_json::to_vec(&vec![feature]).unwrap();
        assert!(json_exact::<FeatureList>(&current).is_ok());
        let text = String::from_utf8(current).unwrap();
        let duplicate = text.replacen("\"kind\":\"x\"", "\"kind\":\"x\",\"kind\":\"x\"", 1);
        assert!(json_exact::<FeatureList>(duplicate.as_bytes()).is_err());
        let unknown = text.replacen("\"kind\":\"x\"", "\"kind\":\"x\",\"unrecognized\":[0]", 1);
        assert!(json_exact::<FeatureList>(unknown.as_bytes()).is_err());
    }
}
