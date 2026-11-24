//! Explicit TA preparation declarations; this is a form adapter, not a kit model.
use collagene_core::{Sequence, analysis::Primer, tailing::TailingBase, topo};
use topo::physical::{self as physical, projection};

#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct Draft {
    pub vector_top: String,
    pub vector_bottom: String,
    pub insert_template: String,
    pub vector_top_length: usize,
    pub vector_bottom_length: usize,
    pub forward_sequence: String,
    pub forward_annealing: usize,
    pub reverse_sequence: String,
    pub reverse_annealing: usize,
    pub preparation_assumption: String,
    pub activation_confirmed: bool,
    pub add_a_tail: bool,
}
fn primer(name: &str, ordered: &str, annealing: usize) -> Result<Primer, String> {
    if ordered.len() > 1024 {
        return Err("An ordered primer exceeds the input text bound.".into());
    }
    let ordered = Sequence::dna(ordered).map_err(|e| e.to_string())?;
    if ordered.is_empty()
        || ordered.len() > collagene_core::primers::MAX_OLIGO
        || annealing == 0
        || annealing > ordered.len()
    {
        return Err(format!(
            "Choose a nonempty 3′ annealing suffix within the ordered {name} primer (at most {} bases).",
            collagene_core::primers::MAX_OLIGO
        ));
    }
    let boundary = ordered.len() - annealing;
    Primer::new(
        name,
        &ordered.as_str()[boundary..],
        &ordered.as_str()[..boundary],
    )
    .map_err(|e| e.to_string())
}

/// Roles map to document IDs. No DNA is copied and no activation is inferred.
/// Returned empty transfers are intentional: orientation inspection must return
/// the exact candidate-specific transfer request before projection/target checks.
pub fn ta_request(draft: &Draft) -> Result<(Vec<(String, String)>, topo::Request), String> {
    for id in [
        &draft.vector_top,
        &draft.vector_bottom,
        &draft.insert_template,
    ] {
        if id.trim().is_empty() || id.len() > 1024 || id.chars().any(char::is_control) {
            return Err("Choose all three source documents.".into());
        }
    }
    if !draft.activation_confirmed {
        return Err("Explicitly declare the prepared vector's 5′ phosphate and 3′ phosphotyrosyl chemistry.".into());
    }
    if !draft.add_a_tail {
        return Err("This TA workflow requires an explicitly declared single A addition to both PCR strand 3′ ends.".into());
    }
    if draft.preparation_assumption.trim().is_empty()
        || draft.preparation_assumption.len() > 4096
        || draft
            .preparation_assumption
            .chars()
            .any(|c| c.is_control() && !matches!(c, '\n' | '\r' | '\t'))
    {
        return Err(
            "Describe the actual vector preparation assumption (up to 4,096 bytes).".into(),
        );
    }
    let length = draft.vector_top_length;
    if length < 2 || length != draft.vector_bottom_length || length > 100_000 {
        return Err("This TA preparation form needs equal-length vector strands of 2–100,000 bases, each with one explicit 3′ overhang. Use the general recipe interface for other geometry.".into());
    }
    let forward = primer("Forward", &draft.forward_sequence, draft.forward_annealing)?;
    let reverse = primer("Reverse", &draft.reverse_sequence, draft.reverse_annealing)?;
    let strands = [
        ("vector_top", "declared-vector-top-3prime"),
        ("vector_bottom", "declared-vector-bottom-3prime"),
    ]
    .into_iter()
    .map(|(role, activation)| topo::PreparedStrand {
        id: role.into(),
        segments: vec![topo::SourceInterval {
            source_role: role.into(),
            start: 0,
            length,
            reverse: false,
        }],
        five_prime: physical::FivePrime::Phosphate,
        three_prime: physical::ThreePrime::Phosphotyrosyl {
            enzyme: physical::Enzyme::VacciniaTypeIBTyr274,
            activation_id: activation.into(),
        },
    })
    .collect();
    let request = topo::Request {
        recipe: topo::GenerationRecipe {
            schema_version: 1,
            materials: vec![
                topo::MaterialRecipe::Prepared {
                    declaration_id: "user-declared-ta-vector-1".into(),
                    assumption: draft.preparation_assumption.clone(),
                    strands,
                    pairing: vec![physical::PairSpan {
                        strand_a: 0,
                        start_a: 0,
                        strand_b: 1,
                        start_b: 0,
                        length: length - 1,
                    }],
                },
                topo::MaterialRecipe::ExactPcr {
                    material_id: "explicit-pcr-a-tail-insert".into(),
                    template_role: "insert_template".into(),
                    forward,
                    reverse,
                    tail: Some(TailingBase::A),
                },
            ],
        },
        transfers: vec![],
        engine_limits: physical::Limits::default(),
        projection_options: projection::Options::default(),
    };
    Ok((
        vec![
            ("vector_top".into(), draft.vector_top.clone()),
            ("vector_bottom".into(), draft.vector_bottom.clone()),
            ("insert_template".into(), draft.insert_template.clone()),
        ],
        request,
    ))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn declaration_preserves_ordered_oligos_and_never_selects_a_transfer() {
        let draft = Draft {
            vector_top: "vt".into(),
            vector_bottom: "vb".into(),
            insert_template: "template".into(),
            vector_top_length: 13,
            vector_bottom_length: 13,
            forward_sequence: "GGACGT".into(),
            forward_annealing: 4,
            reverse_sequence: "TTTGCA".into(),
            reverse_annealing: 4,
            preparation_assumption:
                "Independently supplied prepared vector; explicit chemistry assumption".into(),
            activation_confirmed: true,
            add_a_tail: true,
        };
        let (roles, request) = ta_request(&draft).unwrap();
        assert_eq!(roles[0], ("vector_top".into(), "vt".into()));
        assert!(request.transfers.is_empty());
        let topo::MaterialRecipe::ExactPcr {
            forward,
            reverse,
            tail,
            ..
        } = &request.recipe.materials[1]
        else {
            panic!()
        };
        assert_eq!(forward.annealing.as_str(), "ACGT");
        assert_eq!(forward.tail, "GG");
        assert_eq!(reverse.ordered_sequence(), "TTTGCA");
        assert_eq!(*tail, Some(TailingBase::A));
        let mut missing = draft.clone();
        missing.activation_confirmed = false;
        assert!(ta_request(&missing).is_err());
        missing = draft;
        missing.forward_annealing = 7;
        assert!(ta_request(&missing).is_err());
    }
}
