//! One fresh source reconstruction per merge, never one scan per selection.
use super::*;
use collagene_model::{
    Alphabet, Topology,
    annotations::{MappingSegment, project_features_with_topology},
};
use std::collections::BTreeSet;
pub(super) struct Prepared {
    pub id: String,
    pub feature: Feature,
    pub contribution: Contribution,
}
pub(super) fn prepare(
    target: &Molecule,
    request: &Request,
) -> Result<(Vec<Prepared>, bool, usize), Error> {
    let selected = request
        .selections
        .iter()
        .map(|s| s.proposal_id.as_str())
        .collect::<BTreeSet<_>>();
    let mut output = vec![];
    let mut source_work = 0usize;
    let complete = match &request.source {
        Source::Library { library, report } => {
            if report.comparisons > request.limits.max_comparisons {
                return Err(Error::TooLarge);
            }
            let fresh = library::scan_with_comparisons(
                target,
                library,
                &report.options,
                request.limits.max_comparisons,
            )?;
            source_work = fresh.comparisons;
            if &fresh != report {
                return Err(bad("Library merge source differs from fresh scan"));
            }
            for hit in &fresh.hits {
                if !selected.contains(hit.id.as_str()) {
                    continue;
                }
                let ledger = Ledger::from_feature(&hit.feature)?;
                let contribution=ledger.contributions.into_iter().find(|c|matches!(&c.source,SourcePin::Library{proposal_id,..} if proposal_id==&hit.id)).ok_or_else(||bad("Fresh library proposal has no canonical contribution"))?;
                output.push(Prepared {
                    id: hit.id.clone(),
                    feature: hit.feature.clone(),
                    contribution,
                });
            }
            fresh.complete
        }
        Source::Dataset { report } => {
            report.validate(target)?;
            for proposal in &report.proposals {
                if !selected.contains(proposal.id.as_str()) {
                    continue;
                }
                let source = SourcePin::Dataset {
                    source_blake3: report.source_blake3.clone(),
                    format: report.format,
                    seqid: report.seqid.clone(),
                    proposal_id: proposal.id.clone(),
                    source_rows: proposal.source_rows.clone(),
                };
                let contribution = Contribution::new(source, &proposal.feature)?;
                output.push(Prepared {
                    id: proposal.id.clone(),
                    feature: proposal.feature.clone(),
                    contribution,
                });
            }
            true
        }
        Source::Donor {
            molecule,
            event,
            transform,
            indices,
        } => {
            let n = molecule.sequence.len();
            if molecule.sequence.alphabet() != Alphabet::Dna
                || target.sequence.alphabet() != Alphabet::Dna
                || n != target.sequence.len()
                || molecule.topology != target.topology
            {
                return Err(Error::Unsupported("Direct annotation transfer requires equal-length DNA and identical declared topology".into()));
            }
            if transform.origin >= n
                || (molecule.topology == Topology::Linear && transform.origin != 0)
            {
                return Err(bad(
                    "Donor origin requires a valid declared circular transform",
                ));
            }
            // Check the declared full-record correspondence directly, without a
            // matcher or uniqueness inference, before projecting annotations once.
            if n > request.limits.max_comparisons {
                return Err(Error::TooLarge);
            }
            source_work = n;
            let a = molecule.sequence.as_str().as_bytes();
            let b = target.sequence.as_str().as_bytes();
            for (i, expected) in b.iter().enumerate() {
                let oriented = (i + transform.origin) % n;
                let actual = if transform.reverse {
                    collagene_model::complement(a[n - 1 - oriented], Alphabet::Dna)
                } else {
                    a[oriented]
                };
                if actual != *expected {
                    return Err(bad(
                        "Declared donor transform does not match actual target DNA",
                    ));
                }
            }
            let k = transform.origin;
            let map = if transform.reverse {
                let mut map = vec![MappingSegment {
                    source_start: 0,
                    source_end: n - k,
                    destination_start: 0,
                    reverse: true,
                }];
                if k > 0 {
                    map.push(MappingSegment {
                        source_start: n - k,
                        source_end: n,
                        destination_start: n - k,
                        reverse: true,
                    });
                }
                map
            } else {
                let mut map = vec![MappingSegment {
                    source_start: k,
                    source_end: n,
                    destination_start: 0,
                    reverse: false,
                }];
                if k > 0 {
                    map.push(MappingSegment {
                        source_start: 0,
                        source_end: k,
                        destination_start: n - k,
                        reverse: false,
                    });
                }
                map
            };
            let projection = project_features_with_topology(molecule, &map, n, target.topology)?;
            for &index in indices {
                let id = format!("donor-{index}");
                if !selected.contains(id.as_str()) {
                    continue;
                }
                let original = molecule
                    .features
                    .get(index)
                    .ok_or_else(|| bad("Selected donor annotation is unavailable"))?;
                collagene_model::annotations::biological_path::BiologicalPath::from_location(
                    &original.location,
                    n,
                    collagene_model::annotations::biological_path::PathLimits::default(),
                )?;
                // One explicit complete source feature only; no min/max recovery
                // of uncertain, dropped or multiply emitted annotations.
                let mut found = projection
                    .source_feature_indices
                    .iter()
                    .enumerate()
                    .filter(|(_, source)| **source == index);
                let (destination, _) = found.next().ok_or_else(|| {
                    Error::Unsupported("Selected donor feature is not exactly projectable".into())
                })?;
                if found.next().is_some() {
                    return Err(Error::Unsupported(
                        "Selected donor feature has multiple projection occurrences".into(),
                    ));
                }
                let mut feature = projection.features[destination].clone();
                // Prior donor contributions remain in the exact source event;
                // copying them here would imply independently admitted sources.
                feature.qualifiers.retain(|(key, _)| key != QUALIFIER);
                let source = SourcePin::Donor {
                    event: event.clone(),
                    feature_index: index,
                    annotation_revision: crate::construct_qc::annotation_revision(original)?,
                    transform: transform.clone(),
                };
                let contribution = Contribution::new(source, &feature)?;
                output.push(Prepared {
                    id,
                    feature,
                    contribution,
                });
            }
            true
        }
    };
    if output.len() != request.selections.len() {
        return Err(bad(
            "Selected merge proposal is unavailable in the exact source",
        ));
    }
    if source_work > request.limits.max_comparisons {
        return Err(Error::TooLarge);
    }
    Ok((output, complete, source_work))
}
