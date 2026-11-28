//! Literal sequence-risk candidates, never measured modification annotations.
use crate::protein_liabilities as native;
use collagene_model::{Error, Feature, Location, Molecule};
pub fn candidate(
    source: &Molecule,
    options: &native::Options,
    hit: &native::Hit,
) -> Result<Feature, Error> {
    native::validate_hit(source, options, hit).map_err(|e| Error::Format(e.to_string()))?;
    let context =
        native::request_revision(source, options).map_err(|e| Error::Format(e.to_string()))?;
    let certainty = match hit.certainty {
        native::Certainty::CertainPattern => "certain sequence-pattern match",
        native::Certainty::PossiblePattern => "possible sequence-pattern match",
    };
    Ok(Feature {
        kind: "misc_feature".into(),
        location: Location::range(hit.start, hit.end),
        qualifiers: vec![
            ("label".into(), Some(hit.rule.label().into())),
            (
                "note".into(),
                Some(format!(
                    "Candidate only: {certainty}; not a measured PTM, chemical-state or activity claim. Supplied CDR context, if present, is documentary and not validated numbering."
                )),
            ),
            (
                "inference".into(),
                Some(format!("profile:{}; context:{}", native::ENGINE, context)),
            ),
        ],
        needs_review: true,
    })
}
