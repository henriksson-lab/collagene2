//! Terminal processing previews and freshly derived product documents.
use super::{
    ENGINE_VERSION,
    pcr::project_pcr_annotations,
    projection::{
        ProjectedDocument, ReferenceProjection, project_into, region_map, reverse_reference_map,
    },
};
use crate::{
    chemistry::{self, Duplex},
    sequence::analysis::{self, Primer},
};
use collagene_model::annotations::MappingSegment;
use collagene_model::{Alphabet, Error, Molecule, Topology};
use serde::{Deserialize, Serialize};
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct RepairPreview {
    pub source_revision: String,
    pub enzymes: Vec<String>,
    pub fragment_index: usize,
    pub method: chemistry::EndRepairMethod,
    pub molecule: Molecule,
    pub repair: chemistry::EndRepairReport,
    pub findings: Vec<String>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AnnealPreview {
    pub top_revision: String,
    pub bottom_revision: String,
    pub annealing: chemistry::AnnealReport,
}
fn bounded_physical_preview(value: &impl Serialize) -> Result<(), Error> {
    struct Count(usize);
    impl std::io::Write for Count {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.0 = self
                .0
                .checked_add(bytes.len())
                .filter(|n| *n <= chemistry::MAX_CHEMISTRY_REPORT_BYTES)
                .ok_or_else(|| {
                    std::io::Error::other("physical preview exceeds normal message budget")
                })?;
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    serde_json::to_writer(&mut Count(0), value).map_err(|_| Error::TooLarge)
}
pub(super) fn carry_unmapped_evidence(source: &Molecule, product: &mut Molecule) {
    for (key, value) in &source.metadata {
        if key.starts_with("annotation_mapping:unmapped") {
            product
                .metadata
                .entry(key.clone())
                .or_insert_with(|| value.clone());
        }
    }
}
pub(crate) fn repaired_document(
    source: &Molecule,
    report: &chemistry::EndRepairReport,
    name: String,
    findings: &mut Vec<String>,
) -> Result<ProjectedDocument, Error> {
    let mut product = report.duplex.molecule(name)?;
    let map = report
        .mappings
        .iter()
        .filter(|m| {
            m.source_strand == chemistry::Strand::Top
                && m.product_strand == chemistry::Strand::Top
                && !m.synthesized
        })
        .map(|m| MappingSegment {
            source_start: m.source_start,
            source_end: m.source_end,
            destination_start: m.product_start,
            reverse: m.reverse,
        })
        .collect::<Vec<_>>();
    let native = project_into(source, &map, &mut product, findings)?;
    carry_unmapped_evidence(source, &mut product);
    product
        .metadata
        .insert("end_repair:source_revision".into(), source.revision());
    product
        .metadata
        .insert("end_repair:method".into(), format!("{:?}", report.method));
    product
        .metadata
        .insert("end_repair:engine_version".into(), ENGINE_VERSION.into());
    product.validate()?;
    let introduced = report
        .mappings
        .iter()
        .filter(|m| m.product_strand == chemistry::Strand::Top && m.synthesized)
        .map(|m| (m.product_start, m.product_end))
        .collect();
    Ok(ProjectedDocument {
        molecule: product,
        projections: vec![ReferenceProjection {
            map,
            native,
            introduced,
        }],
    })
}
/// Standalone preview reconstructs the digest rather than accepting an imported
/// duplex summary. The returned molecule is reference/evidence, not physical authority.
pub fn simulate_end_repair(
    source: &Molecule,
    enzymes: &[String],
    fragment_index: usize,
    method: chemistry::EndRepairMethod,
) -> Result<RepairPreview, Error> {
    source.validate()?;
    if enzymes.is_empty() || enzymes.len() > 128 {
        return Err(Error::Format(
            "end-repair preview requires 1–128 explicit restriction enzymes".into(),
        ));
    }
    let catalog = analysis::starter_enzymes();
    let selected = enzymes
        .iter()
        .map(|name| {
            catalog
                .iter()
                .find(|enzyme| enzyme.name.eq_ignore_ascii_case(name))
                .cloned()
                .ok_or_else(|| {
                    Error::Unsupported(format!("enzyme {name} is not in the embedded catalog"))
                })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let initial = Duplex::from_molecule(source)?;
    let fragments = chemistry::digest(&initial, &selected)?;
    let fragment = fragments.get(fragment_index).ok_or_else(|| {
        Error::Format(format!(
            "digest has {} fragments; selected fragment index {fragment_index} is unavailable",
            fragments.len()
        ))
    })?;
    let mut findings = Vec::new();
    let mut intermediate = fragment.duplex.molecule("Reconstructed digest fragment")?;
    project_into(
        source,
        &region_map(
            source.sequence.len(),
            fragment.top_start,
            fragment.duplex.top.len(),
            0,
        )?,
        &mut intermediate,
        &mut findings,
    )?;
    carry_unmapped_evidence(source, &mut intermediate);
    let repair = chemistry::repair_ends(&fragment.duplex, method)?;
    let mut molecule = repaired_document(
        &intermediate,
        &repair,
        format!("{}_repaired_{}", source.name, fragment_index + 1),
        &mut findings,
    )?
    .molecule;
    molecule.metadata.insert(
        "end_repair:original_source_revision".into(),
        source.revision(),
    );
    molecule.metadata.insert(
        "end_repair:enzymes".into(),
        serde_json::to_string(enzymes).map_err(|e| Error::Format(e.to_string()))?,
    );
    molecule.metadata.insert(
        "end_repair:fragment_index".into(),
        fragment_index.to_string(),
    );
    findings.extend(repair.warnings.iter().cloned());
    crate::protocol::save_document(&molecule)?;
    let preview = RepairPreview {
        source_revision: source.revision(),
        enzymes: enzymes.to_vec(),
        fragment_index,
        method,
        molecule,
        repair,
        findings,
    };
    bounded_physical_preview(&preview)?;
    Ok(preview)
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub enum TailSubstrate {
    Digest {
        enzymes: Vec<String>,
        fragment_index: usize,
        blunt_before: bool,
    },
    ExactPcr {
        forward: Primer,
        reverse: Primer,
        template_strand: chemistry::Strand,
    },
}
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TailPreview {
    pub source_revision: String,
    pub recipe: TailSubstrate,
    pub base: crate::cloning::ligation::tailing::TailingBase,
    pub molecule: Molecule,
    pub tailing: crate::cloning::ligation::tailing::TailingReport,
    pub findings: Vec<String>,
}
pub(super) fn tailed_document(
    source: &Molecule,
    report: &crate::cloning::ligation::tailing::TailingReport,
    name: String,
    findings: &mut Vec<String>,
) -> Result<ProjectedDocument, Error> {
    let mut product = report.duplex.molecule(name)?;
    let maps = report
        .mappings
        .iter()
        .filter(|m| {
            m.source_strand == chemistry::Strand::Top
                && m.product_strand == chemistry::Strand::Top
                && !m.synthesized
        })
        .map(|m| MappingSegment {
            source_start: m.source_start,
            source_end: m.source_end,
            destination_start: m.product_start,
            reverse: m.reverse,
        })
        .collect::<Vec<_>>();
    let native = project_into(source, &maps, &mut product, findings)?;
    carry_unmapped_evidence(source, &mut product);
    product
        .metadata
        .insert("tailing:source_revision".into(), source.revision());
    product
        .metadata
        .insert("tailing:base".into(), report.base.as_str().into());
    product
        .metadata
        .insert("tailing:engine_version".into(), ENGINE_VERSION.into());
    product.validate()?;
    let introduced = report
        .additions
        .iter()
        .filter(|addition| addition.strand == chemistry::Strand::Top)
        .map(|addition| (addition.product_start, addition.product_end))
        .collect();
    Ok(ProjectedDocument {
        molecule: product,
        projections: vec![ReferenceProjection {
            map: maps,
            native,
            introduced,
        }],
    })
}
/// Readonly evidence reconstructed from an explicit digest or unique exact PCR.
/// The returned molecule metadata never establishes physical authority in a later plan.
pub fn simulate_tailing(
    source: &Molecule,
    recipe: &TailSubstrate,
    base: crate::cloning::ligation::tailing::TailingBase,
) -> Result<TailPreview, Error> {
    source.validate()?;
    let mut findings = Vec::new();
    let (intermediate, duplex) = match recipe {
        TailSubstrate::Digest {
            enzymes,
            fragment_index,
            blunt_before,
        } => {
            if enzymes.is_empty() || enzymes.len() > 128 {
                return Err(Error::Format(
                    "tailing preview needs1–128 explicit digest enzymes".into(),
                ));
            }
            let catalog = analysis::starter_enzymes();
            let selected = enzymes
                .iter()
                .map(|name| {
                    catalog
                        .iter()
                        .find(|e| e.name.eq_ignore_ascii_case(name))
                        .cloned()
                        .ok_or_else(|| {
                            Error::Unsupported(format!(
                                "enzyme {name} is absent from the embedded catalog"
                            ))
                        })
                })
                .collect::<Result<Vec<_>, _>>()?;
            let initial = Duplex::from_molecule(source)?;
            let fragments = chemistry::digest(&initial, &selected)?;
            let fragment = fragments
                .get(*fragment_index)
                .ok_or_else(|| Error::Format("selected digest fragment is unavailable".into()))?;
            let mut molecule = fragment.duplex.molecule("Tailing digest substrate")?;
            project_into(
                source,
                &region_map(
                    source.sequence.len(),
                    fragment.top_start,
                    fragment.duplex.top.len(),
                    0,
                )?,
                &mut molecule,
                &mut findings,
            )?;
            carry_unmapped_evidence(source, &mut molecule);
            if *blunt_before {
                let report =
                    chemistry::repair_ends(&fragment.duplex, chemistry::EndRepairMethod::Blunt)?;
                let repaired = repaired_document(
                    &molecule,
                    &report,
                    "Blunted tailing substrate".into(),
                    &mut findings,
                )?
                .molecule;
                (repaired, report.duplex)
            } else {
                (molecule, fragment.duplex.clone())
            }
        }
        TailSubstrate::ExactPcr {
            forward,
            reverse,
            template_strand,
        } => {
            if *template_strand != chemistry::Strand::Top {
                return Err(Error::Unsupported("Original tail-preview PCR uses explicit Top; saved Bottom needs native recipe replay".into()));
            }
            if chemistry::has_unmapped_physical_evidence(source) {
                return Err(Error::Unsupported("PCR tailing source has imported physical evidence; reconstruct its originating graph".into()));
            }
            let matches = analysis::simulate_pcr_geometry(source, forward, reverse)?;
            if matches.len() != 1 {
                return Err(Error::Unsupported(format!(
                    "PCR tailing requires exactly one exact product; found {}",
                    matches.len()
                )));
            }
            let geometry = &matches[0];
            let duplex = Duplex::exact_pcr(&geometry.product.sequence)?;
            let mut product = duplex.molecule("Exact PCR tailing substrate")?;
            project_pcr_annotations(source, geometry, &mut product, &mut findings)?;
            (product, duplex)
        }
    };
    let tailing = crate::cloning::ligation::tailing::tail_blunt_duplex(&duplex, base)?;
    let mut molecule = tailed_document(
        &intermediate,
        &tailing,
        format!("{}_{}-tailed", source.name, base.as_str()),
        &mut findings,
    )?
    .molecule;
    molecule
        .metadata
        .insert("tailing:original_source_revision".into(), source.revision());
    molecule.metadata.insert(
        "tailing:recipe".into(),
        serde_json::to_string(recipe).map_err(|e| Error::Format(e.to_string()))?,
    );
    findings.extend(tailing.warnings.iter().cloned());
    let preview = TailPreview {
        source_revision: source.revision(),
        recipe: recipe.clone(),
        base,
        molecule,
        tailing,
        findings,
    };
    bounded_physical_preview(&preview)?;
    Ok(preview)
}
pub(crate) fn declared_oligo(source: &Molecule) -> Result<(), Error> {
    crate::formats::fidelity::declared_oligo_guard(source)?;
    source.validate()?;
    if source.sequence.alphabet() != Alphabet::Dna || source.topology != Topology::Linear {
        return Err(Error::Unsupported(
            "annealing requires declared linear DNA oligo inputs".into(),
        ));
    }
    if source.metadata.contains_key("duplex_model")
        || source
            .metadata
            .get("source:legacy:chemistry_unmapped")
            .is_some_and(|value| value == "true")
        || source
            .metadata
            .get("source:snapgene:chemistry_unmapped")
            .is_some_and(|value| value == "true")
    {
        return Err(Error::Unsupported("imported physical intermediates cannot be flattened into ssDNA oligos; use independently declared oligo inputs".into()));
    }
    if source
        .metadata
        .get("genbank.molecule_type")
        .is_some_and(|value| {
            let value = value.to_ascii_uppercase();
            value.starts_with("DS-") || value.starts_with("MS-")
        })
        || source
            .original_genbank
            .as_ref()
            .and_then(|text| text.lines().find(|line| line.starts_with("LOCUS")))
            .is_some_and(|line| {
                line.split_ascii_whitespace().any(|word| {
                    let word = word.to_ascii_uppercase();
                    word.starts_with("DS-") || word.starts_with("MS-")
                })
            })
        || source.metadata.contains_key("source:snapgene:dna_flags")
    {
        return Err(Error::Unsupported(
            "an explicitly double/mixed-stranded input is not a declared ssDNA oligo".into(),
        ));
    }
    Ok(())
}
pub fn simulate_oligo_annealing(
    top: &Molecule,
    bottom: &Molecule,
    options: &chemistry::AnnealOptions,
) -> Result<AnnealPreview, Error> {
    declared_oligo(top)?;
    declared_oligo(bottom)?;
    let annealing = chemistry::anneal_oligos(&top.sequence, &bottom.sequence, options)?;
    let preview = AnnealPreview {
        top_revision: top.revision(),
        bottom_revision: bottom.revision(),
        annealing,
    };
    bounded_physical_preview(&preview)?;
    Ok(preview)
}
pub(crate) fn annealed_document_inner(
    top: &Molecule,
    bottom: &Molecule,
    report: &chemistry::AnnealReport,
    candidate_index: usize,
    findings: &mut Vec<String>,
) -> Result<ProjectedDocument, Error> {
    let candidate = report
        .candidates
        .get(candidate_index)
        .ok_or_else(|| Error::Format("selected annealing candidate is unavailable".into()))?;
    let mut product = candidate.duplex.molecule("Annealed oligos")?;
    let top_map = vec![MappingSegment {
        source_start: 0,
        source_end: top.sequence.len(),
        destination_start: 0,
        reverse: false,
    }];
    let top_native = project_into(top, &top_map, &mut product, findings)?;
    let end = candidate.duplex.bottom_start + bottom.sequence.len() as i64;
    let bottom_map = reverse_reference_map(
        bottom.sequence.len(),
        end,
        0,
        product.sequence.len(),
        Topology::Linear,
    )?;
    let bottom_native = project_into(bottom, &bottom_map, &mut product, findings)?;
    carry_unmapped_evidence(top, &mut product);
    carry_unmapped_evidence(bottom, &mut product);
    product
        .metadata
        .insert("anneal:top_revision".into(), top.revision());
    product
        .metadata
        .insert("anneal:bottom_revision".into(), bottom.revision());
    product
        .metadata
        .insert("anneal:candidate_id".into(), candidate.id.clone());
    product
        .metadata
        .insert("anneal:engine_version".into(), ENGINE_VERSION.into());
    product.validate()?;
    crate::protocol::save_document(&product)?;
    Ok(ProjectedDocument {
        molecule: product,
        projections: vec![
            ReferenceProjection {
                map: top_map,
                native: top_native,
                introduced: vec![],
            },
            ReferenceProjection {
                map: bottom_map,
                native: bottom_native,
                introduced: vec![],
            },
        ],
    })
}
/// Reconstruct standalone read-only preview evidence; no imported report can
/// provide physical authority or edited annotation coordinates.
pub fn annealed_document(
    top: &Molecule,
    bottom: &Molecule,
    preview: &AnnealPreview,
    candidate_index: usize,
) -> Result<Molecule, Error> {
    if preview.top_revision != top.revision() || preview.bottom_revision != bottom.revision() {
        return Err(Error::Format(
            "annealing preview source revision changed".into(),
        ));
    }
    let rebuilt = simulate_oligo_annealing(top, bottom, &preview.annealing.options)?;
    bounded_physical_preview(preview)?;
    if rebuilt != *preview {
        return Err(Error::Format(
            "annealing preview differs from independently reconstructed evidence".into(),
        ));
    }
    annealed_document_inner(
        top,
        bottom,
        &rebuilt.annealing,
        candidate_index,
        &mut Vec::new(),
    )
    .map(|projected| projected.molecule)
}

#[cfg(test)]
mod receipt_tests {
    use super::*;
    use collagene_model::{Feature, Location};
    fn annotate(source: &mut Molecule, label: &str) {
        source.features.push(Feature {
            kind: "misc_feature".into(),
            location: Location::range(0, source.sequence.len()),
            qualifiers: vec![("label".into(), Some(label.into()))],
            needs_review: false,
        });
    }
    #[test]
    fn repair_and_tail_receipts_bind_actual_retained_maps_and_top_synthesis() {
        // A single cut of this circular substrate creates TWO known cut
        // termini. A linear upstream fragment would retain its unknown outer
        // end, and repairing the cut end must not establish that outer3OH.
        let source = Molecule::dna(
            "Literal circular EcoRI substrate",
            "AAGAATTCTT",
            Topology::Circular,
        )
        .unwrap();
        let fragment = chemistry::digest(
            &Duplex::from_molecule(&source).unwrap(),
            &[analysis::starter_enzymes()
                .into_iter()
                .find(|e| e.name == "EcoRI")
                .unwrap()],
        )
        .unwrap()
        .remove(0);
        // Motif starts2; EcoRI cuts original top3/bottom7. The linearized
        // top is source[3:] + source[:3], with a four-base5′ stagger.
        assert_eq!(fragment.top_start, 3);
        assert_eq!(fragment.duplex.top.as_str(), "AATTCTTAAG");
        assert_eq!(fragment.duplex.bottom_start, 4);
        assert_eq!(fragment.duplex.left, chemistry::Terminus::CLEAVED);
        assert_eq!(fragment.duplex.right, chemistry::Terminus::CLEAVED);
        let mut intermediate = fragment.duplex.molecule("Actual cut intermediate").unwrap();
        annotate(&mut intermediate, "Retained top");
        let repair =
            chemistry::repair_ends(&fragment.duplex, chemistry::EndRepairMethod::Blunt).unwrap();
        let repaired =
            repaired_document(&intermediate, &repair, "Actual repair".into(), &mut vec![]).unwrap();
        assert_eq!(repaired.projections.len(), 1);
        let synthesized = repair
            .mappings
            .iter()
            .filter(|m| m.product_strand == chemistry::Strand::Top && m.synthesized)
            .map(|m| (m.product_start, m.product_end))
            .collect::<Vec<_>>();
        assert_eq!(repaired.projections[0].introduced, synthesized);
        assert!(
            !synthesized.is_empty(),
            "literal staggered circular EcoRI fragment requires top fill-in"
        );
        assert_eq!(
            repaired.projections[0]
                .native
                .projection
                .source_feature_indices,
            vec![0]
        );
        // The right-hand bottom overhang templates AATT at displayed top
        // positions10..14; repair retains the genuinely known cut-end3OHs.
        assert_eq!(synthesized, [(10, 14)]);
        assert_eq!(repair.duplex.top.as_str(), "AATTCTTAAGAATT");
        assert_eq!(repair.duplex.left, chemistry::Terminus::CLEAVED);
        assert_eq!(repair.duplex.right, chemistry::Terminus::CLEAVED);
        let tail = crate::cloning::ligation::tailing::tail_blunt_duplex(
            &repair.duplex,
            crate::cloning::ligation::tailing::TailingBase::A,
        )
        .unwrap();
        let tailed =
            tailed_document(&repaired.molecule, &tail, "Actual tail".into(), &mut vec![]).unwrap();
        assert_eq!(tailed.projections.len(), 1);
        assert_eq!(
            tailed.projections[0].introduced,
            vec![(repair.duplex.top.len(), repair.duplex.top.len() + 1)]
        );
        for (feature, index) in tailed.projections[0]
            .native
            .projection
            .features
            .iter()
            .zip(&tailed.projections[0].native.destination_feature_indices)
        {
            assert_eq!(feature, &tailed.molecule.features[*index]);
        }
    }
    #[test]
    fn annealing_receipts_keep_top_then_bottom_source_ordinals_from_native_projection() {
        let mut top = Molecule::dna("Top oligo", "ACGTTGCA", Topology::Linear).unwrap();
        let mut bottom = Molecule::dna("Bottom oligo", "TGCAACGT", Topology::Linear).unwrap();
        annotate(&mut top, "Top annotation");
        annotate(&mut bottom, "Bottom annotation");
        let mut options = chemistry::AnnealOptions::default();
        options.min_paired_bases = 8;
        let preview = simulate_oligo_annealing(&top, &bottom, &options).unwrap();
        let index = preview
            .annealing
            .candidates
            .iter()
            .position(|c| c.duplex.top.as_str() == "ACGTTGCA" && c.duplex.bottom_start == 0)
            .unwrap();
        let product =
            annealed_document_inner(&top, &bottom, &preview.annealing, index, &mut vec![]).unwrap();
        assert_eq!(product.projections.len(), 2);
        assert!(!product.projections[0].map[0].reverse);
        assert!(product.projections[1].map[0].reverse);
        for receipt in &product.projections {
            assert!(receipt.introduced.is_empty());
            assert_eq!(receipt.native.projection.source_feature_indices, vec![0]);
        }
        let top_original = product.projections[0].native.destination_feature_indices[0];
        let top_final = product.projections[1].native.previous_feature_indices[top_original];
        assert_eq!(
            product.molecule.features[top_final].label(),
            "Top annotation"
        );
        let bottom_final = product.projections[1].native.destination_feature_indices[0];
        assert_eq!(
            product.molecule.features[bottom_final].label(),
            "Bottom annotation"
        );
    }
}
