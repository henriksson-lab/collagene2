//! Conflict review and actual atomic application share one prepared source execution.
use super::*;
use collagene_model::annotations::biological_path::{BiologicalPath, PathLimits};
struct Work {
    used: usize,
    limit: usize,
    relations: usize,
    relation_limit: usize,
}
impl Work {
    fn charge(&mut self) -> bool {
        if self.used == self.limit {
            false
        } else {
            self.used += 1;
            true
        }
    }
}
fn relation(
    a: &Feature,
    b: &Feature,
    n: usize,
    work: &mut Work,
) -> Result<Option<RelationKind>, Error> {
    if !work.charge() {
        return Err(Error::TooLarge);
    }
    let paths = (
        BiologicalPath::from_location(&a.location, n, PathLimits::default()),
        BiologicalPath::from_location(&b.location, n, PathLimits::default()),
    );
    let (a_path, b_path) = match paths {
        (Ok(a), Ok(b)) => (a, b),
        _ => return Ok(Some(RelationKind::Unresolved)),
    };
    let mut equal = a_path.len() == b_path.len();
    let mut opposite = equal;
    let mut overlap = false;
    for i in 0..a_path.len() {
        if !work.charge() {
            return Err(Error::TooLarge);
        }
        let call = a_path.call_at(i).unwrap();
        equal &= b_path.call_at(i) == Some(call);
        if opposite {
            let expected = b_path.call_at(b_path.len() - 1 - i).unwrap();
            opposite = expected.position == call.position && expected.reverse != call.reverse;
        }
        overlap |= b_path.ordinal_of(call).is_some()
            || b_path
                .ordinal_of(
                    collagene_model::annotations::biological_path::BiologicalCall {
                        position: call.position,
                        reverse: !call.reverse,
                    },
                )
                .is_some();
    }
    if equal {
        if a == b {
            Ok(Some(RelationKind::ExactDuplicate))
        } else if provenance::same_meaning(a, b) {
            Ok(Some(RelationKind::EquivalentAnnotation))
        } else {
            Ok(Some(RelationKind::SameGeometryConflict))
        }
    } else if opposite {
        Ok(Some(RelationKind::OppositeOrientation))
    } else if overlap {
        Ok(Some(RelationKind::Overlap))
    } else {
        Ok(None)
    }
}
fn review(target: &Molecule, request: &Request) -> Result<Preview, Error> {
    let context_revision = revision(target, request)?;
    let (prepared, source_complete, source_work) = sources::prepare(target, request)?;
    let mut work = Work {
        used: source_work,
        limit: request.limits.max_comparisons,
        relations: 0,
        relation_limit: request.limits.max_relations,
    };
    let mut candidates = vec![];
    for selection in &request.selections {
        let proposal = prepared
            .iter()
            .find(|p| p.id == selection.proposal_id)
            .ok_or_else(|| bad("Fresh merge proposal missing"))?;
        let mut relations = vec![];
        let mut complete = true;
        for (index, feature) in target.features.iter().enumerate() {
            let kind = match relation(&proposal.feature, feature, target.sequence.len(), &mut work)
            {
                Ok(kind) => kind,
                Err(Error::TooLarge) => {
                    complete = false;
                    break;
                }
                Err(error) => return Err(error),
            };
            if let Some(kind) = kind {
                if work.relations == work.relation_limit {
                    complete = false;
                    break;
                }
                work.relations += 1;
                relations.push(Relation {
                    target_index: index,
                    annotation_revision: crate::construct_qc::annotation_revision(feature)?,
                    kind,
                });
            }
        }
        candidates.push(Candidate {
            proposal_id: proposal.id.clone(),
            feature: proposal.feature.clone(),
            contribution: proposal.contribution.clone(),
            decision: selection.decision.clone(),
            relations,
            relations_complete: complete,
        });
    }
    let complete = source_complete && candidates.iter().all(|c| c.relations_complete);
    let mut output = Preview {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        target_revision: target.revision(),
        source_complete,
        candidates,
        comparisons: work.used,
        complete,
        report_revision: String::new(),
        warnings: vec![
            "Annotation source contributions and declared donor transforms are documentary evidence; function, uniqueness, sample identity and event ownership are not authenticated by this molecule-only assessment. All appended annotations require review. Overlap or equal labels never establish annotation identity.".into(),
        ],
    };
    if !complete {
        output.warnings.push("Source search or target relationship review is incomplete; omitted relationships are not absent.".into());
    }
    output.report_revision = output.identity()?;
    output.validate_request(target, request)?;
    Ok(output)
}
pub(super) fn preview(target: &Molecule, request: &Request) -> Result<Preview, Error> {
    review(target, request)
}
pub(super) fn apply(
    target: &Molecule,
    request: &Request,
) -> Result<(Applied, crate::sequence::editing::AnnotationIdentity), Error> {
    let review = review(target, request)?;
    // Validate every explicit decision and compute prospective attribution on
    // borrowed target ordinals before publishing any changed molecule.
    let mut changes: Vec<(usize, Feature)> = vec![];
    let mut additions: Vec<Feature> = vec![];
    let mut actions = vec![];
    let mut extra_work = Work {
        used: review.comparisons,
        limit: request.limits.max_comparisons,
        relations: 0,
        relation_limit: request.limits.max_relations,
    };
    for candidate in &review.candidates {
        match &candidate.decision {
            Decision::Skip => actions.push(ActionReceipt {
                proposal_id: candidate.proposal_id.clone(),
                decision: Decision::Skip,
                destination_index: None,
                contribution_added: false,
            }),
            Decision::Append => {
                if !candidate.relations_complete {
                    return Err(Error::Unsupported(
                        "Incomplete target relationship review cannot authorize append".into(),
                    ));
                }
                if candidate.relations.iter().any(|r| {
                    matches!(
                        r.kind,
                        RelationKind::ExactDuplicate | RelationKind::EquivalentAnnotation
                    )
                }) {
                    return Err(bad(
                        "Equivalent annotation already exists; choose explicit Skip or AddProvenance",
                    ));
                }
                for prior in &additions {
                    if matches!(
                        relation(
                            &candidate.feature,
                            prior,
                            target.sequence.len(),
                            &mut extra_work
                        )?,
                        Some(RelationKind::ExactDuplicate | RelationKind::EquivalentAnnotation)
                    ) {
                        return Err(bad("Two Append decisions describe equivalent annotations"));
                    }
                }
                let mut feature = candidate.feature.clone();
                feature.needs_review = true;
                Ledger::add_to(&mut feature, &candidate.contribution)?;
                let destination = target.features.len() + additions.len();
                additions.push(feature);
                actions.push(ActionReceipt {
                    proposal_id: candidate.proposal_id.clone(),
                    decision: candidate.decision.clone(),
                    destination_index: Some(destination),
                    contribution_added: true,
                });
            }
            Decision::AddProvenance {
                target_index,
                annotation_revision,
            } => {
                if !candidate.relations_complete {
                    return Err(Error::Unsupported(
                        "Incomplete target relationship review cannot authorize attribution".into(),
                    ));
                }
                let existing = target
                    .features
                    .get(*target_index)
                    .ok_or_else(|| bad("Attribution target ordinal is unavailable"))?;
                if &crate::construct_qc::annotation_revision(existing)? != annotation_revision {
                    return Err(bad("Attribution target annotation pin changed"));
                }
                if !candidate.relations.iter().any(|r| {
                    r.target_index == *target_index
                        && matches!(
                            r.kind,
                            RelationKind::ExactDuplicate | RelationKind::EquivalentAnnotation
                        )
                }) {
                    return Err(bad(
                        "Attribution requires exact ordered geometry and compatible existing meaning",
                    ));
                }
                bounded_json::preflight(existing, request.limits.max_output_bytes)?;
                let position = changes.iter().position(|(index, _)| index == target_index);
                let mut feature = position
                    .map(|i| changes.remove(i).1)
                    .unwrap_or_else(|| existing.clone());
                let added = Ledger::add_to(&mut feature, &candidate.contribution)?;
                changes.push((*target_index, feature));
                actions.push(ActionReceipt {
                    proposal_id: candidate.proposal_id.clone(),
                    decision: candidate.decision.clone(),
                    destination_index: Some(*target_index),
                    contribution_added: added,
                });
            }
        }
    }
    // Borrowed aggregate preflight includes every planned feature before cloning
    // the target. The final frame/JSON guard includes actual shared raw evidence.
    let mut feature_budget = collagene_model::bounds::ModelBudget::annotation_edit();
    for feature in &additions {
        feature_budget.features(std::slice::from_ref(feature))?;
    }
    for (_, feature) in &changes {
        feature_budget.features(std::slice::from_ref(feature))?;
    }
    let dataset_original = match &request.source {
        Source::Dataset { report } if actions.iter().any(|a| a.contribution_added) => {
            Some(&report.original_base64)
        }
        _ => None,
    };
    // Conservative complete prospective output: unchanged target plus planned
    // changes, one actual review, actions and at most one shared external file.
    // Duplicating replaced feature text in this borrowed count is conservative.
    bounded_json::preflight(
        &(
            target,
            &additions,
            &changes,
            &review,
            &actions,
            dataset_original,
        ),
        request.limits.max_output_bytes.saturating_sub(1024),
    )?;
    let mut molecule = target.clone();
    let mut identity =
        crate::sequence::editing::AnnotationIdentity::retained(target.features.len());
    for (index, feature) in changes {
        molecule.features[index] = feature;
    }
    for feature in additions {
        identity.added.push(molecule.features.len());
        molecule.features.push(feature);
    }
    if actions.iter().any(|a| a.contribution_added) {
        if let Source::Dataset { report } = &request.source {
            crate::formats::annotations::retain_merge_source(&mut molecule, report)?;
        }
    }
    molecule.validate()?;
    crate::protocol::save_document(&molecule)?;
    let receipt = MergeReceipt {
        context_revision: review.context_revision.clone(),
        actions,
        retained_original_features: target.features.len(),
        added_feature_indices: identity.added.clone(),
        product_revision: molecule.revision(),
        comparisons: extra_work.used,
    };
    let applied = Applied {
        molecule,
        review,
        receipt,
    };
    applied.validate_binding(&applied.review.context_revision)?;
    bounded_json::preflight(&applied, request.limits.max_output_bytes)?;
    Ok((applied, identity))
}
