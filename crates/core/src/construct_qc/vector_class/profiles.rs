//! Installed data definitions for declared research transfer-plasmid graphs.
use super::*;

pub const RAAV_ID: &str = "raav-transfer-declared-graph-2026";
pub const SIN_LENTI_ID: &str = "sin-lentiviral-transfer-declared-graph-2026";

data! {
    pub struct Source {
        pub id: String,
        pub url: String,
        pub locator: String,
        pub revision: String,
        pub reviewed_on: String,
        pub rights: String,
    }
}
data! {
    #[derive(Copy)]
    pub enum Orientation {
        SameAsAxis,
        OppositeToAxis,
        Either,
    }
}
data! {
    pub struct RoleRule {
        pub id: String,
        pub description: String,
        /// Counts qualified nodes in the explicitly declared graph. This is not
        /// exhaustive biological copy number or proof that an unbound role is absent.
        pub minimum: usize,
        pub maximum: usize,
        pub orientation: Orientation,
        pub within_boundaries: bool,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, String, 8>")]
        pub source_ids: Vec<String>,
    }
}
data! {
    pub struct OrderRule {
        pub upstream_role: String,
        pub downstream_role: String,
        pub minimum_gap: usize,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub maximum_gap: Option<usize>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, String, 8>")]
        pub source_ids: Vec<String>,
    }
}
data! {
    pub struct NumericContext {
        pub id: String,
        pub metric: String,
        pub minimum: usize,
        pub maximum: usize,
        pub unit: String,
        pub population: String,
        pub interpretation: String,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, String, 8>")]
        pub source_ids: Vec<String>,
    }
}
data! {
    pub struct Definition {
        pub id: String,
        pub edition: String,
        pub name: String,
        pub scope: String,
        pub start_role: String,
        pub end_role: String,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, RoleRule, 16>")]
        pub roles: Vec<RoleRule>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, OrderRule, 32>")]
        pub order: Vec<OrderRule>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, NumericContext, 8>")]
        pub numeric_context: Vec<NumericContext>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Source, 8>")]
        pub sources: Vec<Source>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, String, 16>")]
        pub limitations: Vec<String>,
    }
}
data! {
    pub struct ProfileRef {
        pub id: String,
        pub identity: String,
    }
}
data! {
    pub struct Profile {
        pub reference: ProfileRef,
        pub definition: Definition,
    }
}
impl Definition {
    fn validate(&self) -> Result<(), Error> {
        bounded(self, 64 * 1024)?;
        if !text(&self.id, 128)
            || !text(&self.edition, 128)
            || !text(&self.name, 256)
            || !text(&self.scope, 2048)
            || self.roles.is_empty()
            || self.roles.len() > 16
            || self.order.len() > 32
            || self.numeric_context.len() > 8
            || self.sources.is_empty()
            || self.sources.len() > 8
            || self.limitations.len() > 16
            || self.limitations.iter().any(|line| !text(line, 2048))
        {
            return Err(Error::TooLarge);
        }
        let mut sources = BTreeSet::new();
        for source in &self.sources {
            if !text(&source.id, 128)
                || !sources.insert(source.id.as_str())
                || !source.url.starts_with("https://")
                || !text(&source.url, 2048)
                || !text(&source.locator, 1024)
                || !text(&source.revision, 512)
                || !text(&source.reviewed_on, 32)
                || !text(&source.rights, 1024)
            {
                return Err(bad("invalid installed source provenance"));
            }
        }
        let evidence = |ids: &[String]| {
            !ids.is_empty() && ids.len() <= 8 && ids.iter().all(|id| sources.contains(id.as_str()))
        };
        let mut roles = BTreeSet::new();
        for role in &self.roles {
            if !text(&role.id, 128)
                || !roles.insert(role.id.as_str())
                || !text(&role.description, 1024)
                || role.maximum > 64
                || role.minimum > role.maximum
                || !evidence(&role.source_ids)
            {
                return Err(bad("invalid installed role declaration"));
            }
        }
        if !roles.contains(self.start_role.as_str()) || !roles.contains(self.end_role.as_str()) {
            return Err(bad("installed boundary roles are absent"));
        }
        for order in &self.order {
            if !roles.contains(order.upstream_role.as_str())
                || !roles.contains(order.downstream_role.as_str())
                || order.minimum_gap > MAX_BASES
                || order
                    .maximum_gap
                    .is_some_and(|maximum| maximum < order.minimum_gap || maximum > MAX_BASES)
                || !evidence(&order.source_ids)
            {
                return Err(bad("invalid installed directed order rule"));
            }
        }
        for context in &self.numeric_context {
            if !text(&context.id, 128)
                || !text(&context.metric, 256)
                || context.minimum > context.maximum
                || context.maximum > MAX_BASES
                || !text(&context.unit, 32)
                || !text(&context.population, 2048)
                || !text(&context.interpretation, 2048)
                || !evidence(&context.source_ids)
            {
                return Err(bad("invalid descriptive numeric context"));
            }
        }
        Ok(())
    }
}
impl Profile {
    pub fn validate(&self) -> Result<(), Error> {
        self.definition.validate()?;
        if resolve(&self.reference)? != *self {
            return Err(bad("profile differs from exact installed factual curation"));
        }
        Ok(())
    }
}
pub fn catalogue() -> Result<Vec<Profile>, Error> {
    let bytes = include_bytes!("data/profiles.json");
    if bytes.len() > 128 * 1024 {
        return Err(Error::TooLarge);
    }
    let definitions: Vec<Definition> = serde_json::from_slice(bytes)
        .map_err(|error| Error::Format(format!("installed vector-class data: {error}")))?;
    if definitions.len() > 16 {
        return Err(Error::TooLarge);
    }
    let mut ids = BTreeSet::new();
    definitions
        .into_iter()
        .map(|definition| {
            definition.validate()?;
            if !ids.insert(definition.id.clone()) {
                return Err(bad("duplicate installed vector profile"));
            }
            let reference = ProfileRef {
                id: definition.id.clone(),
                identity: hash(&(ENGINE_VERSION, &definition))?,
            };
            Ok(Profile {
                reference,
                definition,
            })
        })
        .collect()
}
pub fn catalogue_revision() -> Result<String, Error> {
    hash(&(ENGINE_VERSION, catalogue()?))
}
pub fn resolve(reference: &ProfileRef) -> Result<Profile, Error> {
    catalogue()?
        .into_iter()
        .find(|profile| profile.reference == *reference)
        .ok_or_else(|| bad("unknown or altered installed vector-class profile pin"))
}
