//! PCR preview products and their annotation projection.
use super::{
    ENGINE_VERSION,
    projection::{ReferenceProjection, project_into, region_map},
};
use crate::{
    chemistry::{self, tracked_duplex::TrackedDuplex},
    sequence::analysis::{self, Primer},
};
use collagene_model::{Error, Molecule, Topology};
use serde::{Deserialize, Serialize};
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PcrDocument {
    pub molecule: Molecule,
    pub geometry: analysis::PcrProduct,
    pub findings: Vec<String>,
}

/// Enumerate sequence-reference previews. These are not target verification reports or
/// evidence of physical end chemistry; the immutable source remains the parent document.
pub fn simulate_pcr_documents(
    source: &Molecule,
    forward: &Primer,
    reverse: &Primer,
    template_strand: chemistry::Strand,
) -> Result<Vec<PcrDocument>, Error> {
    build_pcr_documents(
        reference_frame(source, template_strand)?,
        forward,
        reverse,
        None,
        None,
    )
    .map(|(documents, _)| documents)
}
/// Retain each actual native projection by the same geometry index. None means
/// receipt retention exceeded its own bound; scientific products are unaffected.
pub(crate) fn simulate_pcr_documents_with_projections(
    source: &Molecule,
    forward: &Primer,
    reverse: &Primer,
    template_strand: chemistry::Strand,
    budget: Option<&mut analysis::ExactSearchBudget>,
) -> Result<(Vec<PcrDocument>, Option<Vec<ReferenceProjection>>), Error> {
    build_pcr_documents(
        reference_frame(source, template_strand)?,
        forward,
        reverse,
        Some(RECEIPT_BYTES),
        budget,
    )
}
fn reference_frame(
    source: &Molecule,
    strand: chemistry::Strand,
) -> Result<crate::cloning::pcr::source::Template<'_>, Error> {
    if strand != chemistry::Strand::Top {
        return Err(Error::Unsupported("Original reference PCR uses explicit Top; Bottom requires freshly reconstructed native source strands".into()));
    }
    crate::cloning::pcr::source::Template::reference(source)
}
/// Only actual fresh native state supplies a selected Bottom template frame.
pub(crate) fn simulate_pcr_documents_from_fresh(
    source: &Molecule,
    forward: &Primer,
    reverse: &Primer,
    tracked: &TrackedDuplex,
    template_strand: chemistry::Strand,
) -> Result<Vec<PcrDocument>, Error> {
    build_pcr_documents(
        crate::cloning::pcr::source::Template::fresh(source, tracked, template_strand)?,
        forward,
        reverse,
        None,
        None,
    )
    .map(|(documents, _)| documents)
}
pub(crate) fn simulate_pcr_documents_from_fresh_with_projections(
    source: &Molecule,
    forward: &Primer,
    reverse: &Primer,
    tracked: &TrackedDuplex,
    template_strand: chemistry::Strand,
    budget: Option<&mut analysis::ExactSearchBudget>,
) -> Result<(Vec<PcrDocument>, Option<Vec<ReferenceProjection>>), Error> {
    build_pcr_documents(
        crate::cloning::pcr::source::Template::fresh(source, tracked, template_strand)?,
        forward,
        reverse,
        Some(RECEIPT_BYTES),
        budget,
    )
}
/// Share the enclosing native replay's actual matching/materialization ledger.
pub(crate) fn simulate_pcr_documents_with_budget(
    source: &Molecule,
    forward: &Primer,
    reverse: &Primer,
    template_strand: chemistry::Strand,
    budget: &mut analysis::ExactSearchBudget,
) -> Result<Vec<PcrDocument>, Error> {
    build_pcr_documents(
        reference_frame(source, template_strand)?,
        forward,
        reverse,
        None,
        Some(budget),
    )
    .map(|(documents, _)| documents)
}
pub(crate) fn simulate_pcr_documents_from_fresh_with_budget(
    source: &Molecule,
    forward: &Primer,
    reverse: &Primer,
    tracked: &TrackedDuplex,
    template_strand: chemistry::Strand,
    budget: &mut analysis::ExactSearchBudget,
) -> Result<Vec<PcrDocument>, Error> {
    build_pcr_documents(
        crate::cloning::pcr::source::Template::fresh(source, tracked, template_strand)?,
        forward,
        reverse,
        None,
        Some(budget),
    )
    .map(|(documents, _)| documents)
}
fn build_pcr_documents(
    frame: crate::cloning::pcr::source::Template<'_>,
    forward: &Primer,
    reverse: &Primer,
    capture_limit: Option<usize>,
    budget: Option<&mut analysis::ExactSearchBudget>,
) -> Result<(Vec<PcrDocument>, Option<Vec<ReferenceProjection>>), Error> {
    let source = frame.reference;
    let geometries = match budget {
        Some(budget) => frame.geometries_with_budget(forward, reverse, budget)?.0,
        None => frame.geometries(forward, reverse)?,
    };
    let mut total = 0usize;
    let mut documents = Vec::new();
    let mut projections = capture_limit.map(|_| Vec::new());
    let mut receipt_bytes = 0usize;
    for (index, geometry) in geometries.into_iter().enumerate() {
        let mut molecule = Molecule::dna(
            format!("{}_PCR_{}", source.name, index + 1),
            geometry.product.sequence.as_str(),
            Topology::Linear,
        )?;
        let mut findings = vec![format!(
            "PCR coordinates use the actual {:?} source frame; complementary-strand coverage was checked without filling, repairing or cropping the template.",
            frame.strand
        )];
        let projection =
            project_pcr_frame_annotations(frame, &geometry, &mut molecule, &mut findings)?;
        molecule
            .metadata
            .insert("pcr:template_revision".into(), source.revision());
        molecule
            .metadata
            .insert("pcr:source_frame".into(), frame.context_revision()?);
        molecule
            .metadata
            .insert("pcr:template_strand".into(), format!("{:?}", frame.strand));
        molecule.metadata.insert(
            "pcr:primers".into(),
            serde_json::to_string(&(forward, reverse)).map_err(|e| Error::Format(e.to_string()))?,
        );
        molecule.metadata.insert("pcr:model".into(), "Exact annealing reference preview; alternative geometries remain separate. No target verification or physical chemistry is asserted.".into());
        molecule
            .metadata
            .insert("pcr:engine_version".into(), ENGINE_VERSION.into());
        let size = serde_json::to_vec(&molecule)
            .map_err(|e| Error::Format(e.to_string()))?
            .len();
        total = total.checked_add(size).ok_or(Error::TooLarge)?;
        if total > crate::formats::sequence::MAX_IMPORT_BYTES {
            return Err(Error::TooLarge);
        }
        documents.push(PcrDocument {
            molecule,
            geometry: geometry.product,
            findings,
        });
        if let Some(limit) = capture_limit.filter(|_| projections.is_some()) {
            match receipt_size(&projection, limit.saturating_sub(receipt_bytes)) {
                Some(bytes) => {
                    receipt_bytes += bytes;
                    if let Some(projections) = &mut projections {
                        projections.push(projection);
                    }
                }
                None => projections = None,
            }
        }
    }
    Ok((documents, projections))
}
const RECEIPT_BYTES: usize = 4 * 1024 * 1024;
fn receipt_size(receipt: &ReferenceProjection, maximum: usize) -> Option<usize> {
    struct Count {
        used: usize,
        maximum: usize,
    }
    impl std::io::Write for Count {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.used = self
                .used
                .checked_add(bytes.len())
                .filter(|n| *n <= self.maximum)
                .ok_or_else(|| std::io::Error::other("PCR receipt retention bound"))?;
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let projection = &receipt.native.projection;
    let mut count = Count { used: 0, maximum };
    serde_json::to_writer(
        &mut count,
        &(
            &receipt.map,
            &receipt.introduced,
            &receipt.native.destination_feature_indices,
            &receipt.native.previous_feature_indices,
            &projection.features,
            &projection.source_feature_indices,
            &projection.warnings,
            &projection.unmapped_features,
            &projection.unmapped_source_feature_indices,
            projection.dropped_features,
            &projection.dropped_source_feature_indices,
        ),
    )
    .ok()?;
    Some(count.used)
}

pub(crate) fn project_pcr_annotations(
    source: &Molecule,
    geometry: &analysis::PcrGeometryProduct,
    product: &mut Molecule,
    findings: &mut Vec<String>,
) -> Result<ReferenceProjection, Error> {
    if geometry.product.template_strand != chemistry::Strand::Top {
        return Err(Error::Format(
            "Bottom PCR annotation projection requires the exact native frame".into(),
        ));
    }
    let map = region_map(
        source.sequence.len(),
        geometry.product.forward_start,
        geometry.template_length,
        geometry.forward_tail_length,
    )?;
    let native = project_into(source, &map, product, findings)?;
    Ok(ReferenceProjection {
        map,
        native,
        introduced: tail_ranges(
            product.sequence.len(),
            geometry.forward_tail_length,
            geometry.reverse_tail_length,
        )?,
    })
}

pub(crate) fn project_pcr_frame_annotations(
    frame: crate::cloning::pcr::source::Template<'_>,
    geometry: &analysis::PcrGeometryProduct,
    product: &mut Molecule,
    findings: &mut Vec<String>,
) -> Result<ReferenceProjection, Error> {
    let map = frame.annotation_map(geometry)?;
    if frame.strand == chemistry::Strand::Bottom {
        findings.push("Bottom-frame PCR maps original top annotations only where real source calls correspond; bottom-only extensions have no fabricated top coordinates.".into());
    }
    let native = project_into(frame.reference, &map, product, findings)?;
    Ok(ReferenceProjection {
        map,
        native,
        introduced: tail_ranges(
            product.sequence.len(),
            geometry.forward_tail_length,
            geometry.reverse_tail_length,
        )?,
    })
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MutagenesisPreview {
    pub report: crate::cloning::pcr::mutagenesis::PcrReport,
    pub documents: Vec<Molecule>,
    pub findings: Vec<String>,
}
pub fn simulate_mutagenesis_preview(
    source: &Molecule,
    forward: &Primer,
    reverse: &Primer,
    options: &crate::cloning::pcr::mutagenesis::PcrOptions,
) -> Result<MutagenesisPreview, Error> {
    if chemistry::has_unmapped_physical_evidence(source) {
        return Err(Error::Unsupported(
            "PCR preview of physical-intermediate evidence requires its originating graph".into(),
        ));
    }
    let report = crate::cloning::pcr::mutagenesis::simulate_pcr(source, forward, reverse, options)?;
    project_mutagenesis_preview(source, forward, reverse, options, report)
}
/// Project the one freshly computed native report, preserving its full search evidence.
pub(crate) fn project_mutagenesis_preview(
    source: &Molecule,
    forward: &Primer,
    reverse: &Primer,
    options: &crate::cloning::pcr::mutagenesis::PcrOptions,
    report: crate::cloning::pcr::mutagenesis::PcrReport,
) -> Result<MutagenesisPreview, Error> {
    let mut documents = Vec::new();
    let mut findings = Vec::new();
    let mut bytes = 0usize;
    for (i, geometry) in report.products.iter().enumerate() {
        let mut product = Molecule::dna(
            format!("{}_mutant_PCR_{}", source.name, i + 1),
            geometry.sequence.as_str(),
            Topology::Linear,
        )?;
        project_mutant_pcr(source, geometry, &mut product, &mut findings)?;
        product
            .metadata
            .insert("pcr:template_revision".into(), source.revision());
        product.metadata.insert(
            "pcr:mismatch_options".into(),
            serde_json::to_string(options).map_err(|e| Error::Format(e.to_string()))?,
        );
        product.metadata.insert(
            "pcr:primers".into(),
            serde_json::to_string(&(forward, reverse)).map_err(|e| Error::Format(e.to_string()))?,
        );
        product.metadata.insert("pcr:model".into(),"Intentional-mismatch reference preview; no target verification, amplification or physical end chemistry is asserted. Retain the source document and report.".into());
        product
            .metadata
            .insert("pcr:engine_version".into(), ENGINE_VERSION.into());
        bytes = bytes.saturating_add(
            serde_json::to_vec(&product)
                .map_err(|e| Error::Format(e.to_string()))?
                .len(),
        );
        if bytes > 16 * 1024 * 1024 {
            return Err(Error::TooLarge);
        }
        documents.push(product);
    }
    Ok(MutagenesisPreview {
        report,
        documents,
        findings,
    })
}
pub(crate) fn project_mutant_pcr(
    source: &Molecule,
    geometry: &crate::cloning::pcr::mutagenesis::MutantPcrProduct,
    product: &mut Molecule,
    findings: &mut Vec<String>,
) -> Result<ReferenceProjection, Error> {
    let map = region_map(
        source.sequence.len(),
        geometry.forward_start,
        geometry.template_length,
        geometry.forward_tail_length,
    )?;
    let mut native = project_into(source, &map, product, findings)?;
    let mut changed = geometry
        .substitutions
        .iter()
        .map(|change| change.product_position)
        .collect::<Vec<_>>();
    changed.sort_unstable();
    changed.dedup();
    for feature in &mut product.features {
        if feature.location.spans().iter().any(|(a, b, _)| {
            changed
                .get(changed.partition_point(|position| *position < *a))
                .is_some_and(|position| *position < *b)
        }) {
            feature.needs_review = true;
            let mut evidence = Vec::new();
            feature.qualifiers.retain(|(key, value)| {
                if key == "translation" || key == "codon_start" {
                    evidence.push(format!(
                        "Source assertion /{key}={}",
                        value.as_deref().unwrap_or("")
                    ));
                    false
                } else {
                    true
                }
            });
            feature.qualifiers.push(("note".into(),Some("Primer-introduced changes overlap this annotation; review its sequence and function.".into())));
            for note in evidence {
                feature.qualifiers.push(("note".into(), Some(note)));
            }
        }
    }
    product.validate()?;
    // Annotation review edits are part of the same projection. Bind the actual
    // retained feature, including its removed source translation assertions.
    for (feature, destination) in native
        .projection
        .features
        .iter_mut()
        .zip(&native.destination_feature_indices)
    {
        *feature = product
            .features
            .get(*destination)
            .ok_or_else(|| {
                Error::Format("Mutant PCR projection destination ordinal differs".into())
            })?
            .clone();
    }
    Ok(ReferenceProjection {
        map,
        native,
        introduced: tail_ranges(
            product.sequence.len(),
            geometry.forward_tail_length,
            geometry.reverse_tail_length,
        )?,
    })
}

fn tail_ranges(
    length: usize,
    forward: usize,
    reverse: usize,
) -> Result<Vec<(usize, usize)>, Error> {
    if forward
        .checked_add(reverse)
        .is_none_or(|tails| tails > length)
    {
        return Err(Error::Coordinates(
            "PCR tail intervals exceed the actual product".into(),
        ));
    }
    let mut ranges = Vec::new();
    if forward > 0 {
        ranges.push((0, forward));
    }
    if reverse > 0 {
        ranges.push((length - reverse, length));
    }
    Ok(ranges)
}

#[cfg(test)]
mod receipt_tests {
    use super::*;
    use collagene_model::{Feature, Location};
    fn annotated() -> Molecule {
        let mut source = Molecule::dna(
            "Actual primer frame",
            "ACGTTGCATGATCCGATGCA",
            Topology::Linear,
        )
        .unwrap();
        source.features.push(Feature {
            kind: "CDS".into(),
            location: Location::range(0, 12),
            qualifiers: vec![
                ("label".into(), Some("Source CDS".into())),
                ("translation".into(), Some("source assertion".into())),
            ],
            needs_review: false,
        });
        source
    }
    #[test]
    fn exact_pcr_receipt_uses_executed_map_ordinals_and_only_literal_tails() {
        let source = annotated();
        let forward = Primer::new("forward", "ACGTTGCA", "GG").unwrap();
        let reverse = Primer::new("reverse", "TGCATCGG", "TTT").unwrap();
        let geometry = analysis::simulate_pcr_geometry(&source, &forward, &reverse)
            .unwrap()
            .into_iter()
            .find(|g| g.product.forward_start == 0 && g.template_length == 20)
            .unwrap();
        let frame = crate::cloning::pcr::source::Template::reference(&source).unwrap();
        let mut product = Molecule::dna(
            "Actual PCR",
            geometry.product.sequence.as_str(),
            Topology::Linear,
        )
        .unwrap();
        let receipt =
            project_pcr_frame_annotations(frame, &geometry, &mut product, &mut vec![]).unwrap();
        assert_eq!(receipt.introduced, vec![(0, 2), (22, 25)]);
        assert_eq!(receipt.map.len(), 1);
        assert_eq!(
            (
                receipt.map[0].source_start,
                receipt.map[0].source_end,
                receipt.map[0].destination_start
            ),
            (0, 20, 2)
        );
        assert_eq!(receipt.native.projection.source_feature_indices, vec![0]);
        assert_eq!(receipt.native.destination_feature_indices, vec![0]);
        assert_eq!(receipt.native.projection.features[0], product.features[0]);
    }
    #[test]
    fn receipt_retention_refusal_keeps_exact_native_products_and_geometry_indices() {
        let source = annotated();
        let forward = Primer::new("forward", "ACGTTGCA", "GG").unwrap();
        let reverse = Primer::new("reverse", "TGCATCGG", "TTT").unwrap();
        let (documents, receipts) = simulate_pcr_documents_with_projections(
            &source,
            &forward,
            &reverse,
            chemistry::Strand::Top,
            None,
        )
        .unwrap();
        let receipts = receipts.unwrap();
        assert_eq!(documents.len(), receipts.len());
        for (document, receipt) in documents.iter().zip(&receipts) {
            for (feature, index) in receipt
                .native
                .projection
                .features
                .iter()
                .zip(&receipt.native.destination_feature_indices)
            {
                assert_eq!(feature, &document.molecule.features[*index]);
            }
        }
        let (limited, refused) = build_pcr_documents(
            reference_frame(&source, chemistry::Strand::Top).unwrap(),
            &forward,
            &reverse,
            Some(0),
            None,
        )
        .unwrap();
        assert!(refused.is_none());
        assert_eq!(
            serde_json::to_value(&documents).unwrap(),
            serde_json::to_value(&limited).unwrap()
        );
        let ordinary =
            simulate_pcr_documents(&source, &forward, &reverse, chemistry::Strand::Top).unwrap();
        assert_eq!(
            serde_json::to_value(&documents).unwrap(),
            serde_json::to_value(&ordinary).unwrap()
        );
    }
    #[test]
    fn enclosing_pcr_budget_survives_failed_queries_and_receipt_refusal() {
        let source = annotated();
        let forward = Primer::new("forward", "ACGTTGCA", "GG").unwrap();
        let reverse = Primer::new("reverse", "TGCATCGG", "TTT").unwrap();
        let mut captured_budget = analysis::ExactSearchBudget::new(100_000).unwrap();
        let (captured, receipts) = simulate_pcr_documents_with_projections(
            &source,
            &forward,
            &reverse,
            chemistry::Strand::Top,
            Some(&mut captured_budget),
        )
        .unwrap();
        assert!(receipts.is_some());
        let mut refused_budget = analysis::ExactSearchBudget::new(100_000).unwrap();
        let (refused_products, refused) = build_pcr_documents(
            reference_frame(&source, chemistry::Strand::Top).unwrap(),
            &forward,
            &reverse,
            Some(0),
            Some(&mut refused_budget),
        )
        .unwrap();
        assert!(refused.is_none());
        assert!(captured_budget.used() > 0);
        assert_eq!(captured_budget.used(), refused_budget.used());
        assert_eq!(
            serde_json::to_value(captured).unwrap(),
            serde_json::to_value(refused_products).unwrap()
        );

        let absent = Primer::new("absent", "AAAAAAAAAAAA", "").unwrap();
        let mut shared = analysis::ExactSearchBudget::new(1).unwrap();
        assert!(
            simulate_pcr_documents_with_budget(
                &source,
                &absent,
                &absent,
                chemistry::Strand::Top,
                &mut shared,
            )
            .is_err()
        );
        assert_eq!(shared.remaining(), 0);
        let spent = shared.used();
        assert!(
            simulate_pcr_documents_with_projections(
                &source,
                &forward,
                &reverse,
                chemistry::Strand::Top,
                Some(&mut shared),
            )
            .is_err()
        );
        assert_eq!(shared.used(), spent);
    }
    #[test]
    fn mutant_receipt_retains_actual_post_review_feature_without_source_translation_assertion() {
        let source = annotated();
        let forward = Primer::new("first-base substitution", "CCGTTGCA", "").unwrap();
        let reverse = Primer::new("reverse", "TGCATCGG", "").unwrap();
        let report = crate::cloning::pcr::mutagenesis::simulate_pcr(
            &source,
            &forward,
            &reverse,
            &Default::default(),
        )
        .unwrap();
        let geometry = report
            .products
            .iter()
            .find(|g| g.forward_start == 0 && g.template_length == 20)
            .unwrap();
        assert!(
            geometry
                .substitutions
                .iter()
                .any(|c| c.product_position == 0 && c.product_base == 'C')
        );
        let mut product = Molecule::dna(
            "Actual mismatch PCR",
            geometry.sequence.as_str(),
            Topology::Linear,
        )
        .unwrap();
        let receipt = project_mutant_pcr(&source, geometry, &mut product, &mut vec![]).unwrap();
        assert!(
            receipt.introduced.is_empty(),
            "a substitution is not a new literal tail"
        );
        assert_eq!(receipt.native.projection.source_feature_indices, vec![0]);
        let destination = receipt.native.destination_feature_indices[0];
        let feature = &product.features[destination];
        assert!(feature.needs_review);
        assert!(
            !feature
                .qualifiers
                .iter()
                .any(|(key, _)| key == "translation")
        );
        assert_eq!(&receipt.native.projection.features[0], feature);
        assert!(
            source.features[0]
                .qualifiers
                .iter()
                .any(|(key, _)| key == "translation")
        );
    }
}
