//! Product-scoped host methylase documentation, never molecular modification state.
//! Selecting a record does not assert propagation, complete occupancy, loss of
//! incoming methylation, CpG status, or the absence of other methylase systems.
use crate::{Catalog, Citation, Error, bounded, fail};
use serde::{Deserialize, Serialize};

pub const MAX_HOSTS: usize = 128;
pub const MAX_HOST_BYTES: usize = 32 * 1024;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(rename_all = "snake_case")]
pub enum MethylaseExpectation {
    Active,
    Deficient,
    Unknown,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct HostPreset {
    pub id: String,
    pub name: String,
    pub species: String,
    pub vendor: String,
    pub product_code: String,
    pub host_scope: String,
    pub dam: MethylaseExpectation,
    pub dcm: MethylaseExpectation,
    pub cpg: MethylaseExpectation,
    #[serde(deserialize_with = "bounded_vec::<_, String, 16>")]
    pub qualifiers: Vec<String>,
    #[serde(deserialize_with = "bounded_vec::<_, Citation, 16>")]
    pub citations: Vec<Citation>,
}
impl HostPreset {
    /// Structural bounds only; current catalogue authority requires validate_host.
    pub fn validate(&self) -> Result<(), Error> {
        if !bounded(&self.id, 96)
            || !bounded(&self.name, 128)
            || !bounded(&self.species, 128)
            || !bounded(&self.vendor, 128)
            || !bounded(&self.product_code, 96)
            || !bounded(&self.host_scope, 1024)
            || self.qualifiers.is_empty()
            || self.qualifiers.len() > 16
            || self.qualifiers.iter().any(|s| !bounded(s, 1024))
            || self.citations.is_empty()
            || self.citations.len() > 16
            || self
                .citations
                .iter()
                .any(|c| !bounded(&c.source, 96) || !bounded(&c.locator, 512))
        {
            return fail("invalid documentary source-host record");
        }
        serde_json::to_writer(crate::ByteBudget(MAX_HOST_BYTES), self)
            .map_err(|_| Error("source-host record exceeds byte limit".into()))
    }
}
impl Catalog {
    pub fn hosts(&self) -> &[HostPreset] {
        &self.edition.hosts
    }
    /// Exact catalogue ID lookup. Strain families and derivatives are not aliases.
    pub fn host(&self, id: &str) -> Result<&HostPreset, Error> {
        if !bounded(id, 96) {
            return fail("invalid source-host identity");
        }
        self.hosts()
            .iter()
            .find(|h| h.id == id)
            .ok_or_else(|| Error("source-host identity is not installed in this catalogue".into()))
    }
    /// Value membership includes scope, qualifications and every source citation.
    /// A matching ID/name/hash supplied by a caller cannot substitute new facts.
    pub fn validate_host(&self, host: &HostPreset) -> Result<(), Error> {
        host.validate()?;
        if self.host(&host.id)? != host {
            return fail("source-host record differs from current catalogue");
        }
        Ok(())
    }
}
pub(crate) fn deserialize_hosts<'de, D: serde::Deserializer<'de>>(
    d: D,
) -> Result<Vec<HostPreset>, D::Error> {
    bounded_vec::<D, HostPreset, MAX_HOSTS>(d)
}
fn bounded_vec<'de, D, T, const MAX: usize>(d: D) -> Result<Vec<T>, D::Error>
where
    D: serde::Deserializer<'de>,
    T: Deserialize<'de>,
{
    struct Rows<T, const N: usize>(std::marker::PhantomData<T>);
    impl<'de, T: Deserialize<'de>, const N: usize> serde::de::Visitor<'de> for Rows<T, N> {
        type Value = Vec<T>;
        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            write!(f, "at most {N} records")
        }
        fn visit_seq<A: serde::de::SeqAccess<'de>>(
            self,
            mut seq: A,
        ) -> Result<Self::Value, A::Error> {
            use serde::de::Error;
            if seq.size_hint().is_some_and(|n| n > N) {
                return Err(A::Error::custom("too many host records"));
            }
            let mut rows = Vec::with_capacity(seq.size_hint().unwrap_or(0).min(N));
            for _ in 0..N {
                match seq.next_element()? {
                    Some(row) => rows.push(row),
                    None => return Ok(rows),
                }
            }
            if seq.next_element::<serde::de::IgnoredAny>()?.is_some() {
                return Err(A::Error::custom("too many host records"));
            }
            Ok(rows)
        }
    }
    d.deserialize_seq(Rows::<T, MAX>(std::marker::PhantomData))
}
