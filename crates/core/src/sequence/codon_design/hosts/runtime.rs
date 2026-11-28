//! Compact pinned runtime catalogue. Raw CDS/GFF derivation is maintainer/worker work.
use super::{HostProfile, bad, digest};
use collagene_model::Error;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Reference {
    pub id: String,
    pub identity: String,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Summary {
    pub reference: Reference,
    pub organism: String,
    pub assembly: String,
    pub genetic_code: u8,
    pub admitted_loci: u64,
    pub excluded_loci: u64,
    pub internal_codons: u64,
    pub policy: String,
}
fn mg1655() -> Result<HostProfile, Error> {
    let profile: HostProfile =
        serde_json::from_str(include_str!("data/mg1655/profile.json")).map_err(bad)?;
    profile.validate()?;
    Ok(profile)
}
impl Summary {
    fn from_profile(profile: &HostProfile) -> Result<Self, Error> {
        let excluded_loci = profile
            .observed_loci
            .checked_sub(profile.admitted_loci)
            .ok_or_else(|| bad("invalid admitted/observed locus counts"))?;
        Ok(Self {
            reference: Reference {
                id: profile.dataset.id.clone(),
                identity: profile.identity.clone(),
            },
            organism: profile.dataset.organism.clone(),
            assembly: profile.dataset.assembly.clone(),
            genetic_code: profile.dataset.genetic_code,
            admitted_loci: profile.admitted_loci,
            excluded_loci,
            internal_codons: profile.internal_codons,
            policy: profile.policy.clone(),
        })
    }
}
pub fn catalogue() -> Result<Vec<Summary>, Error> {
    Ok(vec![Summary::from_profile(&mg1655()?)?])
}
pub fn resolve(reference: &Reference) -> Result<HostProfile, Error> {
    if reference.id.len() > 128 || !digest(&reference.identity) {
        return Err(bad("invalid host profile reference"));
    }
    let profile = mg1655()?;
    if reference.id != profile.dataset.id || reference.identity != profile.identity {
        return Err(Error::Unsupported(
            "unknown current host codon profile identity".into(),
        ));
    }
    Ok(profile)
}
