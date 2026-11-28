//! Shared source-coordinate projection into freshly reconstructed products.
use crate::chemistry;
use collagene_model::annotations::{self, MappingSegment};
use collagene_model::{Error, Molecule, Topology};
pub(super) fn region_map(
    length: usize,
    start: usize,
    count: usize,
    destination: usize,
) -> Result<Vec<MappingSegment>, Error> {
    if start >= length || count > length || count == 0 {
        return Err(Error::Coordinates(
            "annotation source interval is outside one source turn".into(),
        ));
    }
    let first = count.min(length - start);
    let mut map = vec![MappingSegment {
        source_start: start,
        source_end: start + first,
        destination_start: destination,
        reverse: false,
    }];
    if first < count {
        map.push(MappingSegment {
            source_start: 0,
            source_end: count - first,
            destination_start: destination + first,
            reverse: false,
        });
    }
    Ok(map)
}
/// Old top base x is on the reversed bottom strand at boundary E-x.
/// Once ligated, neighboring top segments can restore its pairing, so clip
/// against the final reference rather than the original paired intersection.
pub(super) fn reverse_reference_map(
    source_length: usize,
    bottom_end: i64,
    destination: usize,
    product_length: usize,
    topology: Topology,
) -> Result<Vec<MappingSegment>, Error> {
    let end = destination as i64 + bottom_end;
    let n = product_length as i64;
    if topology == Topology::Linear {
        let start = (end - n).max(0).min(source_length as i64) as usize;
        let stop = end.max(0).min(source_length as i64) as usize;
        return Ok(if start < stop {
            vec![MappingSegment {
                source_start: start,
                source_end: stop,
                destination_start: (end - stop as i64) as usize,
                reverse: true,
            }]
        } else {
            vec![]
        });
    }
    if source_length > product_length {
        return Err(Error::Unsupported(
            "reversed strand exceeds one final circular reference turn".into(),
        ));
    }
    let mut map = Vec::new();
    let mut start = 0;
    while start < source_length {
        let boundary = (end - start as i64).rem_euclid(n) as usize;
        let upper = if boundary == 0 {
            product_length
        } else {
            boundary
        };
        let count = (source_length - start).min(upper);
        map.push(MappingSegment {
            source_start: start,
            source_end: start + count,
            destination_start: upper - count,
            reverse: true,
        });
        start += count;
    }
    Ok(map)
}

/// Native annotation projection plus the actual deduplication ordinal map.
/// Existing callers may discard this receipt; fresh ancestry capture consumes it.
pub(crate) struct NativeProjection {
    pub projection: annotations::Projection,
    pub destination_feature_indices: Vec<usize>,
    pub previous_feature_indices: Vec<usize>,
}
/// Coordinates and annotations emitted together by one actual native helper.
pub(crate) struct ReferenceProjection {
    pub map: Vec<MappingSegment>,
    pub native: NativeProjection,
    pub introduced: Vec<(usize, usize)>,
}
/// An actual native source operand ordinal, never a name/content match.
pub(crate) struct SourceProjection {
    pub input: usize,
    pub reference: ReferenceProjection,
}
#[derive(Clone, Debug)]
pub(crate) struct TopoScope {
    pub occurrence: u32,
    pub binding: std::sync::Arc<super::ancestry::OriginBinding>,
}
pub(crate) enum GraphProjectionSource {
    CurrentInput {
        index: usize,
    },
    RetainedGraph {
        role: String,
        scopes: Vec<TopoScope>,
        retained_index: Option<usize>,
        molecule: Molecule,
    },
}
pub(crate) struct GraphProjections {
    pub retained: Vec<crate::material::run::RetainedReplay>,
    pub sources: Vec<GraphProjectionSource>,
    pub projections: Vec<SourceProjection>,
}
pub(crate) struct ProjectedDocument {
    pub molecule: Molecule,
    pub projections: Vec<ReferenceProjection>,
}
pub(crate) fn project_into(
    source: &Molecule,
    map: &[MappingSegment],
    product: &mut Molecule,
    findings: &mut Vec<String>,
) -> Result<NativeProjection, Error> {
    project_into_with_collapse(source, map, product, findings, false)
}
/// Shared source-top coordinate projection for a freshly reconstructed digest fragment.
/// The caller owns cleavage eligibility; this helper creates no physical authority.
pub(crate) fn project_digest_annotations(
    source: &Molecule,
    fragment: &chemistry::DigestFragment,
    product: &mut Molecule,
    findings: &mut Vec<String>,
) -> Result<ReferenceProjection, Error> {
    let map = region_map(
        source.sequence.len(),
        fragment.top_start,
        fragment.duplex.top.len(),
        0,
    )?;
    let native = project_into(source, &map, product, findings)?;
    Ok(ReferenceProjection {
        map,
        native,
        introduced: vec![],
    })
}
pub(super) fn project_into_with_collapse(
    source: &Molecule,
    map: &[MappingSegment],
    product: &mut Molecule,
    findings: &mut Vec<String>,
    allow_homology_collapse: bool,
) -> Result<NativeProjection, Error> {
    if source.features.is_empty() {
        return Ok(NativeProjection {
            projection: annotations::Projection {
                features: vec![],
                source_feature_indices: vec![],
                warnings: vec![],
                unmapped_features: vec![],
                unmapped_source_feature_indices: vec![],
                dropped_features: 0,
                dropped_source_feature_indices: vec![],
            },
            destination_feature_indices: vec![],
            previous_feature_indices: (0..product.features.len()).collect(),
        });
    }
    let projection = if allow_homology_collapse {
        annotations::project_collapsing_features(
            source,
            map,
            product.sequence.len(),
            product.topology,
        )?
    } else {
        annotations::project_features_with_topology(
            source,
            map,
            product.sequence.len(),
            product.topology,
        )?
    };
    let old_length = product.features.len();
    product.features.extend(projection.features.iter().cloned());
    if !projection.unmapped_features.is_empty() {
        product.metadata.insert(
            format!("annotation_mapping:unmapped:{}", source.revision()),
            serde_json::to_string(&projection.unmapped_features)
                .map_err(|e| Error::Format(e.to_string()))?,
        );
    }
    for warning in projection.warnings.iter().take(100) {
        findings.push(warning.chars().take(2048).collect());
    }
    if projection.warnings.len() > 100 {
        findings.push("Additional annotation-projection findings omitted; source evidence and review flags remain in the plan.".into());
    }
    // Native first-wins deduplication records every actual old/new ordinal.
    // Earlier receipts are remapped with this exact retain-loop result rather
    // than assuming that an appended feature still has its original ordinal.
    let mut seen = std::collections::BTreeMap::new();
    let mut remap = Vec::with_capacity(product.features.len());
    let mut next = 0usize;
    product.features.retain(|feature| {
        let hash =
            *blake3::hash(&serde_json::to_vec(feature).expect("serializable feature")).as_bytes();
        if let Some(index) = seen.get(&hash) {
            remap.push(*index);
            false
        } else {
            seen.insert(hash, next);
            remap.push(next);
            next += 1;
            true
        }
    });
    product.validate()?;
    Ok(NativeProjection {
        projection,
        destination_feature_indices: remap[old_length..].to_vec(),
        previous_feature_indices: remap[..old_length].to_vec(),
    })
}
