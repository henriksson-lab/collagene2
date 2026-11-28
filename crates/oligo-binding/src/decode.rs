//! Transport admission, including a report-wide aligned-column allocation budget.
use crate::{Column, Hit, MAX_HITS, MAX_OLIGO_BASES, MAX_RETAINED_COLUMNS};
use serde::{
    Deserializer,
    de::{self, SeqAccess, Visitor},
};
use std::{cell::Cell, fmt};

thread_local! { static REMAINING: Cell<Option<usize>> = const { Cell::new(None) }; }
struct Scope(Option<usize>);
impl Scope {
    fn enter() -> Self {
        Self(REMAINING.with(|remaining| {
            let previous = remaining.get();
            if previous.is_none() {
                remaining.set(Some(MAX_RETAINED_COLUMNS));
            }
            previous
        }))
    }
}
impl Drop for Scope {
    fn drop(&mut self) {
        // A nested scope shares its parent's remaining allowance. Only the outer
        // report owns initialization/reset, including early errors or unwinding.
        if self.0.is_none() {
            REMAINING.with(|remaining| remaining.set(None));
        }
    }
}
pub(crate) fn columns<'de, D: Deserializer<'de>>(deserializer: D) -> Result<Vec<Column>, D::Error> {
    struct Columns;
    impl<'de> Visitor<'de> for Columns {
        type Value = Vec<Column>;
        fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
            f.write_str("bounded aligned columns")
        }
        fn visit_seq<A: SeqAccess<'de>>(self, mut sequence: A) -> Result<Self::Value, A::Error> {
            let allowance = REMAINING
                .with(|r| r.get().unwrap_or(MAX_OLIGO_BASES))
                .min(MAX_OLIGO_BASES);
            if sequence.size_hint().is_some_and(|n| n > allowance) {
                return Err(collagene_model::decode::admission(
                    collagene_model::Error::TooLarge,
                ));
            }
            let mut columns = Vec::new();
            while let Some(column) = sequence.next_element::<Column>()? {
                if columns.len() >= MAX_OLIGO_BASES {
                    return Err(collagene_model::decode::admission(
                        collagene_model::Error::TooLarge,
                    ));
                }
                let admitted = REMAINING.with(|remaining| match remaining.get() {
                    Some(0) => false,
                    Some(n) => {
                        remaining.set(Some(n - 1));
                        true
                    }
                    None => true,
                });
                if !admitted {
                    return Err(collagene_model::decode::admission(
                        collagene_model::Error::TooLarge,
                    ));
                }
                columns.push(column);
            }
            Ok(columns)
        }
    }
    deserializer.deserialize_seq(Columns)
}
pub(crate) fn hits<'de, D: Deserializer<'de>>(deserializer: D) -> Result<Vec<Hit>, D::Error> {
    struct Hits;
    impl<'de> Visitor<'de> for Hits {
        type Value = Vec<Hit>;
        fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
            f.write_str("bounded report hits")
        }
        fn visit_seq<A: SeqAccess<'de>>(self, mut sequence: A) -> Result<Self::Value, A::Error> {
            if sequence.size_hint().is_some_and(|n| n > MAX_HITS) {
                return Err(collagene_model::decode::admission(
                    collagene_model::Error::TooLarge,
                ));
            }
            let _scope = Scope::enter();
            let mut hits = Vec::new();
            while let Some(hit) = sequence.next_element::<Hit>()? {
                if hits.len() >= MAX_HITS {
                    return Err(collagene_model::decode::admission(
                        collagene_model::Error::TooLarge,
                    ));
                }
                hits.push(hit);
            }
            Ok(hits)
        }
    }
    deserializer.deserialize_seq(Hits)
}
fn text<'de, D: Deserializer<'de>, const MAX: usize>(deserializer: D) -> Result<String, D::Error> {
    struct Text<const MAX: usize>;
    impl<'de, const MAX: usize> Visitor<'de> for Text<MAX> {
        type Value = String;
        fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "string of at most {MAX} bytes")
        }
        fn visit_str<E: de::Error>(self, value: &str) -> Result<String, E> {
            if value.len() > MAX {
                return Err(collagene_model::decode::admission(
                    collagene_model::Error::TooLarge,
                ));
            }
            Ok(value.to_owned())
        }
        fn visit_string<E: de::Error>(self, value: String) -> Result<String, E> {
            if value.len() > MAX {
                return Err(collagene_model::decode::admission(
                    collagene_model::Error::TooLarge,
                ));
            }
            Ok(value)
        }
    }
    deserializer.deserialize_string(Text::<MAX>)
}
pub(crate) fn identifier<'de, D: Deserializer<'de>>(d: D) -> Result<String, D::Error> {
    text::<D, 128>(d)
}
pub(crate) fn name<'de, D: Deserializer<'de>>(d: D) -> Result<String, D::Error> {
    text::<D, 256>(d)
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde::de::value::{Error, SeqDeserializer, StrDeserializer};
    #[test]
    fn announced_counts_reject_before_elements_are_decoded() {
        // Unit elements are deliberately invalid hits/columns: count preflight
        // must reject without attempting to deserialize the first element.
        assert!(
            hits(SeqDeserializer::<_, Error>::new(
                std::iter::repeat(()).take(MAX_HITS + 1)
            ))
            .is_err()
        );
        assert!(
            columns(SeqDeserializer::<_, Error>::new(
                std::iter::repeat(()).take(MAX_OLIGO_BASES + 1)
            ))
            .is_err()
        );
    }
    #[test]
    fn aggregate_column_allowance_precedes_hit_allocation_and_scope_restores() {
        assert_eq!(REMAINING.with(Cell::get), None);
        {
            let _scope = Scope::enter();
            REMAINING.with(|r| r.set(Some(1)));
            assert!(
                columns(SeqDeserializer::<_, Error>::new(
                    std::iter::repeat(()).take(2)
                ))
                .is_err()
            );
            assert_eq!(REMAINING.with(Cell::get), Some(1));
            {
                let _nested = Scope::enter();
                assert_eq!(REMAINING.with(Cell::get), Some(1));
            }
            assert_eq!(REMAINING.with(Cell::get), Some(1));
        }
        assert_eq!(REMAINING.with(Cell::get), None);
    }
    #[test]
    fn failed_report_decode_does_not_leak_allowance_into_next_report() {
        assert!(hits(SeqDeserializer::<_, Error>::new(std::iter::once(()))).is_err());
        assert_eq!(REMAINING.with(Cell::get), None);
        assert!(
            hits(SeqDeserializer::<_, Error>::new(std::iter::empty::<()>()))
                .unwrap()
                .is_empty()
        );
        assert_eq!(REMAINING.with(Cell::get), None);
    }
    #[test]
    fn identifiers_reject_oversize_before_borrowed_text_copy() {
        let text = "a".repeat(129);
        assert!(identifier(StrDeserializer::<Error>::new(&text)).is_err());
        assert_eq!(
            identifier(StrDeserializer::<Error>::new("o1")).unwrap(),
            "o1"
        );
    }
}
