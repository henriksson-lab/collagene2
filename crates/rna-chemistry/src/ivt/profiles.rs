//! Small independently curated rule set, never a polymerase activity model.
use super::{ENGINE, ProfileReference, hash};
use crate::Result;
use serde::{Deserialize, Serialize};

pub const T7_ID: &str = "t7-consensus-cleancap-ag-template-2026-11-28";
pub const T7_PROMOTER: &str = "TAATACGACTCACTATA";

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Source {
    pub id: String,
    pub url: String,
    pub locator: String,
    pub reviewed_on: String,
    pub revision: String,
    pub rights: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Profile {
    pub reference: ProfileReference,
    pub promoter_dna: String,
    pub cap_ag_initiator_dna: String,
    pub scope: String,
    #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Source, 8>")]
    pub sources: Vec<Source>,
}

impl Profile {
    /// Authenticity is installed exact membership, not a caller-provided URL/hash.
    pub fn validate(&self) -> Result<()> {
        if self != &current()? {
            return Err(crate::Error::Binding);
        }
        Ok(())
    }
}

pub fn profiles() -> Result<Vec<Profile>> {
    Ok(vec![current()?])
}

pub(super) fn current() -> Result<Profile> {
    let sources = vec![
        Source {
            id: "neb-t7-ag-start".into(),
            url: "https://www.neb.com/en/faqs/do-i-need-to-change-my-promoter-sequence-to-use-cleancap-reagent-ag".into(),
            locator: "T7 promoter and initiating bases +1/+2".into(),
            reviewed_on: "2026-11-28".into(),
            revision: "Live supplier FAQ reviewed 2026-11-28; no immutable revision supplied".into(),
            rights: "Independent factual curation; supplier prose/assets not redistributed".into(),
        },
        Source {
            id: "trilink-n7113-ag".into(),
            url: "https://shop.trilinkbiotech.com/shopping/products/ivt-and-capping-reagents/rna-cap-analogs/cleancap-reagent-ag/".into(),
            locator: "N-7113 description and comparison with ARCA; +1/+2 AG".into(),
            reviewed_on: "2026-11-28".into(),
            revision: "Live product page version 2026.11.25 build 5011".into(),
            rights: "Independent factual curation; supplier prose/assets not redistributed; no technology licence implied".into(),
        },
        Source {
            id: "neb-runoff-ends".into(),
            url: "https://www.neb.com/en-gb/faqs/how-to-prepare-plasmid-templates".into(),
            locator: "ssRNA template downstream linearization and blunt/5-prime versus 3-prime overhang guidance".into(),
            reviewed_on: "2026-11-28".into(),
            revision: "Live supplier FAQ reviewed 2026-11-28; no immutable revision supplied".into(),
            rights: "Independent factual curation; supplier prose/assets not redistributed".into(),
        },
        Source {
            id: "conrad-2020-context".into(),
            url: "https://doi.org/10.1038/s42003-020-01167-x".into(),
            locator: "Results/Figure 2: +1..+3 context and 5-prime heterogeneity; Figure 3: SP6 is distinct".into(),
            reviewed_on: "2026-11-28".into(),
            revision: "Conrad et al., Communications Biology 3:439 (2020), published 2020-08-14".into(),
            rights: "CC BY 4.0, authors2020; independently summarized qualitative scope; no parameter table copied".into(),
        },
    ];
    let promoter_dna = T7_PROMOTER.to_owned();
    let cap_ag_initiator_dna = "AG".to_owned();
    let scope = "Declared T7 consensus-promoter template and original CleanCap AG N-7113 initiator compatibility only; no other polymerase/cap analogue, yield or end-homogeneity model".to_owned();
    let identity = hash(&(
        ENGINE,
        T7_ID,
        &promoter_dna,
        &cap_ag_initiator_dna,
        &scope,
        &sources,
    ))?;
    Ok(Profile {
        reference: ProfileReference {
            id: T7_ID.into(),
            identity,
        },
        promoter_dna,
        cap_ag_initiator_dna,
        scope,
        sources,
    })
}
