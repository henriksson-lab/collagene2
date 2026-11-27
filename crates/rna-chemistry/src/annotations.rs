use crate::*;
use collagene_model::Location;
const REFERENCE: &str = "https://www.insdc.org/submitting-standards/feature-table/";
data! {pub struct FeatureReview {pub feature_index:usize,pub expected_feature_revision:String,pub site_entry:usize}}
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct AnnotationProfile {
    pub code: &'static str,
    pub parent: char,
    pub nucleoside: Option<Nucleoside>,
    pub sugar: Option<Sugar>,
}
pub fn annotation_profiles() -> Vec<AnnotationProfile> {
    use Nucleoside::*;
    [
        ("m6a", 'A', Some(N6Methyladenosine), None),
        ("m5c", 'C', Some(FiveMethylcytidine), None),
        ("m5u", 'U', Some(FiveMethyluridine), None),
        ("p", 'U', Some(Pseudouridine), None),
        ("i", 'A', Some(Inosine), None),
        ("cm", 'C', None, Some(Sugar::TwoPrimeOMethyl)),
        ("gm", 'G', None, Some(Sugar::TwoPrimeOMethyl)),
        ("um", 'U', None, Some(Sugar::TwoPrimeOMethyl)),
        ("fm", 'U', Some(Pseudouridine), Some(Sugar::TwoPrimeOMethyl)),
        (
            "tm",
            'U',
            Some(FiveMethyluridine),
            Some(Sugar::TwoPrimeOMethyl),
        ),
    ]
    .into_iter()
    .map(|(code, parent, nucleoside, sugar)| AnnotationProfile {
        code,
        parent,
        nucleoside,
        sugar,
    })
    .collect()
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AnnotationCandidate {
    pub feature_index: usize,
    pub feature_revision: String,
    pub position: Option<usize>,
    pub raw_code: Option<String>,
    pub suggested_site: Option<Site>,
    pub unsupported_reason: Option<String>,
}
pub fn annotation_candidates(source: &Molecule) -> Result<Vec<AnnotationCandidate>> {
    source_revision(source)?;
    candidates(source)
}
fn candidates(source: &Molecule) -> Result<Vec<AnnotationCandidate>> {
    let profiles = annotation_profiles();
    let bases = source.sequence.as_str().as_bytes();
    let mut out = vec![];
    for evidence in retained_evidence(source) {
        let RetainedEvidence::Feature { index } = evidence else {
            continue;
        };
        let f = &source.features[index];
        let digest = hash(b"collagene-rna-feature-1\0", &(index, f), MAX_SOURCE_BYTES)?;
        let values = f
            .qualifiers
            .iter()
            .filter(|(k, _)| k.eq_ignore_ascii_case("mod_base"))
            .collect::<Vec<_>>();
        let raw_code = if values.len() == 1 {
            values[0].1.clone()
        } else {
            None
        };
        let position = match &f.location {
            Location::Range {
                start,
                end,
                before: false,
                after: false,
            } if start.checked_add(1) == Some(*end) => Some(*start),
            _ => None,
        };
        let code = raw_code
            .as_deref()
            .unwrap_or("")
            .trim()
            .to_ascii_lowercase();
        let profile = profiles.iter().find(|p| p.code == code);
        let reason = if f.needs_review {
            Some("Feature was marked for review after editing")
        } else if !f.kind.eq_ignore_ascii_case("modified_base") {
            Some("Only explicit modified_base features can be reviewed")
        } else if position.is_none() {
            Some(
                "Modification needs one exact stored-strand residue; fuzzy/joined/complement locations are unsupported",
            )
        } else if raw_code.is_none() {
            Some("Exactly one nonempty mod_base qualifier is required")
        } else if profile.is_none() {
            Some("Modification code is outside the current INSDC profile")
        } else if bases.get(position.unwrap()).copied() != Some(profile.unwrap().parent as u8) {
            Some("Modification parent contradicts the stored reference base")
        } else {
            None
        };
        let suggested_site = if reason.is_none() {
            let p = profile.unwrap();
            Some(Site {
                index: position.unwrap(),
                nucleoside: p.nucleoside.clone().unwrap_or(Nucleoside::Canonical),
                sugar: p.sugar.clone().unwrap_or(Sugar::Ribose),
                provenance: Provenance {
                    statement: format!(
                        "Review of imported modified_base feature {} with /mod_base={code}; documentary source constraint",
                        index + 1
                    ),
                    reference: Some(REFERENCE.into()),
                },
            })
        } else {
            None
        };
        out.push(AnnotationCandidate {
            feature_index: index,
            feature_revision: digest,
            position,
            raw_code,
            suggested_site,
            unsupported_reason: reason.map(str::to_owned),
        });
    }
    Ok(out)
}
pub(crate) fn validate_reviews(source: &Molecule, d: &Declaration) -> Result<()> {
    if d.annotation_reviews.len() > 4096 {
        return Err(Error::Limit("annotation reviews"));
    }
    if d.annotation_reviews.is_empty() {
        return Ok(());
    }
    let candidates = candidates(source)?
        .into_iter()
        .map(|c| (c.feature_index, c))
        .collect::<std::collections::BTreeMap<_, _>>();
    let profiles = annotation_profiles();
    let mut seen = BTreeSet::new();
    for r in &d.annotation_reviews {
        if !seen.insert(r.feature_index) {
            return Err(Error::Invalid("duplicate feature review"));
        }
        let c = candidates
            .get(&r.feature_index)
            .ok_or(Error::Invalid("review references no modification feature"))?;
        if c.feature_revision != r.expected_feature_revision {
            return Err(Error::Binding);
        }
        let suggested = c.suggested_site.as_ref().ok_or(Error::Invalid(
            "unsupported or uncertain annotation cannot be accepted",
        ))?;
        let actual = d.sites.get(r.site_entry).ok_or(Error::Invalid(
            "review references no declared residue entry",
        ))?;
        let code = c.raw_code.as_ref().unwrap().trim().to_ascii_lowercase();
        let p = profiles.iter().find(|p| p.code == code).unwrap();
        if actual.index != suggested.index
            || p.nucleoside
                .as_ref()
                .is_some_and(|n| n != &actual.nucleoside)
            || p.sugar.as_ref().is_some_and(|s| s != &actual.sugar)
        {
            return Err(Error::Invalid(
                "reviewed annotation conflicts with the declared residue chemistry",
            ));
        }
    }
    Ok(())
}
