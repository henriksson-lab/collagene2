//! Shared encoding contracts, independent of scientific workflow ownership.
use serde::Deserialize;

pub(crate) mod bounded_json;

/// Admit explicit null while letting serde reject a missing required field.
pub(crate) fn required_option<'de, D: serde::Deserializer<'de>, T: Deserialize<'de>>(
    deserializer: D,
) -> Result<Option<T>, D::Error> {
    Option::deserialize(deserializer)
}
