//! PCR preview products and their annotation projection.
use super::{
    ENGINE_VERSION,
    projection::{project_into, region_map},
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
    if template_strand != chemistry::Strand::Top {
        return Err(Error::Unsupported("Original reference PCR uses explicit Top; Bottom requires freshly reconstructed native source strands".into()));
    }
    build_pcr_documents(
        crate::cloning::pcr::source::Template::reference(source)?,
        forward,
        reverse,
    )
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
    )
}
fn build_pcr_documents(
    frame: crate::cloning::pcr::source::Template<'_>,
    forward: &Primer,
    reverse: &Primer,
) -> Result<Vec<PcrDocument>, Error> {
    let source = frame.reference;
    let geometries = frame.geometries(forward, reverse)?;
    let mut total = 0usize;
    let mut documents = Vec::new();
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
    }
    Ok(documents)
}

pub(crate) fn project_pcr_annotations(
    source: &Molecule,
    geometry: &analysis::PcrGeometryProduct,
    product: &mut Molecule,
    findings: &mut Vec<String>,
) -> Result<(), Error> {
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
    project_into(source, &map, product, findings)
}

pub(crate) fn project_pcr_frame_annotations(
    frame: crate::cloning::pcr::source::Template<'_>,
    geometry: &analysis::PcrGeometryProduct,
    product: &mut Molecule,
    findings: &mut Vec<String>,
) -> Result<(), Error> {
    let map = frame.annotation_map(geometry)?;
    if frame.strand == chemistry::Strand::Bottom {
        findings.push("Bottom-frame PCR maps original top annotations only where real source calls correspond; bottom-only extensions have no fabricated top coordinates.".into());
    }
    project_into(frame.reference, &map, product, findings)
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
) -> Result<(), Error> {
    let map = region_map(
        source.sequence.len(),
        geometry.forward_start,
        geometry.template_length,
        geometry.forward_tail_length,
    )?;
    project_into(source, &map, product, findings)?;
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
    product.validate()
}
