//! Text exchange of the actual reviewed product, using the native format writers.
use collagene_core::formats;
use collagene_model::{Molecule, bounds::model_preflight};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Format {
    GenBank,
    Fasta,
}
impl Format {
    pub(crate) fn name(self) -> &'static str {
        match self {
            Self::GenBank => "GenBank",
            Self::Fasta => "FASTA",
        }
    }
}

#[derive(Debug)]
pub(crate) struct Prepared {
    pub source_revision: String,
    pub format: Format,
    pub filename: String,
    pub text: String,
    pub warnings: Vec<String>,
}
impl Prepared {
    pub(crate) fn require_source(&self, source: &Molecule) -> Result<(), String> {
        model_preflight(&[source]).map_err(|e| e.to_string())?;
        if source.revision() != self.source_revision {
            return Err(
                "The actual reconstructed product changed; prepare its export again.".into(),
            );
        }
        Ok(())
    }
    pub(crate) fn loss_report(&self, review_identity: &str) -> Result<Vec<u8>, String> {
        if review_identity.len() > 16 * 1024 {
            return Err("Review identity exceeds the export report limit.".into());
        }
        serde_json::to_vec_pretty(&serde_json::json!({
            "schema_version": 1,
            "scope": "Actual reconstructed product reference; may be diagnostic. Text export does not install material or accept an intended target.",
            "review_identity": review_identity,
            "product_revision": self.source_revision,
            "format": self.format.name(),
            "filename": self.filename,
            "bytes": self.text.len(),
            "warnings": self.warnings,
        }))
        .map_err(|e| e.to_string())
    }
}

pub(crate) fn prepare(source: &Molecule, format: Format) -> Result<Prepared, String> {
    model_preflight(&[source]).map_err(|e| e.to_string())?;
    source.validate().map_err(|e| e.to_string())?;
    let native_format = match format {
        Format::GenBank => formats::TextFormat::GenBank,
        Format::Fasta => formats::TextFormat::Fasta,
    };
    let mut result =
        formats::export_with_report(source, native_format).map_err(|e| e.to_string())?;
    if result.text.len() > formats::MAX_IMPORT_BYTES {
        return Err(
            "Text export exceeds the supported 32 MiB exchange-file limit. Nothing was downloaded."
                .into(),
        );
    }
    result.warnings.push("This actual-product reference may be diagnostic, including a product that differs from the intended target. Text exchange does not retain the complete plan, exact source-event history or physical graph evidence, and does not restore material authority on reimport. Retain the native workspace and full reconstruction report.".into());
    let stem: String = source
        .name
        .chars()
        .take(80)
        .map(|c| {
            if c.is_ascii_alphanumeric() || matches!(c, '-' | '_') {
                c
            } else {
                '_'
            }
        })
        .collect();
    let stem = if stem.is_empty() { "unnamed" } else { &stem };
    let extension = match format {
        Format::GenBank => "gb",
        Format::Fasta => "fasta",
    };
    Ok(Prepared {
        source_revision: source.revision(),
        format,
        filename: format!("actual-product-{stem}.{extension}"),
        text: result.text,
        warnings: result.warnings,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use collagene_core::chemistry;
    use collagene_model::{Alphabet, Feature, Location, Sequence, Topology};

    fn actual_product() -> Molecule {
        let pcr = chemistry::Duplex::exact_pcr(&Sequence::dna("AACCGGTT").unwrap()).unwrap();
        let phosphorylated = pcr.treat(chemistry::EndTreatment::Phosphorylate).unwrap();
        let mut product = chemistry::ligate(&[phosphorylated], true, false)
            .unwrap()
            .molecule("actual / diagnostic")
            .unwrap();
        product.features.push(Feature {
            kind: "misc_feature".into(),
            location: Location::Complement(Box::new(Location::Join(vec![
                Location::range(6, 8),
                Location::range(0, 3),
            ]))),
            qualifiers: vec![(
                "label".into(),
                Some("retained origin-spanning feature".into()),
            )],
            needs_review: true,
        });
        product
    }

    #[test]
    fn genbank_uses_native_annotation_roundtrip_and_unverified_physical_comments() {
        let product = actual_product();
        let original = product.clone();
        let prepared = prepare(&product, Format::GenBank).unwrap();
        let reread = formats::import_text(&prepared.text)
            .unwrap()
            .molecules
            .remove(0);
        assert_eq!(reread.sequence.as_str(), "AACCGGTT");
        assert_eq!(reread.topology, Topology::Circular);
        assert_eq!(reread.features[0].location, product.features[0].location);
        assert_eq!(reread.features[0].label(), product.features[0].label());
        assert!(reread.features[0].needs_review);
        assert!(chemistry::has_unmapped_physical_evidence(&reread));
        assert!(chemistry::Duplex::from_molecule(&reread).is_err());
        assert!(
            prepared
                .warnings
                .iter()
                .any(|w| w.contains("Duplex end states/nicks"))
        );
        assert_eq!(product, original);
        assert!(!prepared.filename.contains('/'));
    }

    #[test]
    fn fasta_discloses_real_losses_and_loss_report_binds_the_actual_product() {
        let product = actual_product();
        let prepared = prepare(&product, Format::Fasta).unwrap();
        let reread = formats::import_text(&prepared.text)
            .unwrap()
            .molecules
            .remove(0);
        assert_eq!(reread.sequence, product.sequence);
        assert_eq!(reread.topology, Topology::Linear);
        assert!(reread.features.is_empty());
        assert!(
            prepared
                .warnings
                .iter()
                .any(|w| w.contains("FASTA exports sequence only"))
        );
        let json: serde_json::Value =
            serde_json::from_slice(&prepared.loss_report("plan+target").unwrap()).unwrap();
        assert_eq!(json["product_revision"], product.revision());
        assert_eq!(json["bytes"], prepared.text.len());
        assert_eq!(
            json["warnings"].as_array().unwrap().len(),
            prepared.warnings.len()
        );
        assert!(json.get("text").is_none());
        let mut changed = product.clone();
        changed.features[0].needs_review = false;
        assert!(prepared.require_source(&changed).is_err());
        assert!(prepared.require_source(&product).is_ok());
    }

    #[test]
    fn native_writer_limitations_and_model_bounds_remain_explicit() {
        let mut protein = Molecule::dna("protein", "AAA", Topology::Linear).unwrap();
        protein.sequence = Sequence::new(Alphabet::Protein, "MKL").unwrap();
        assert!(
            prepare(&protein, Format::GenBank)
                .unwrap_err()
                .contains("GenPept")
        );
        let fasta = prepare(&protein, Format::Fasta).unwrap();
        let reread = formats::import_text(&fasta.text)
            .unwrap()
            .molecules
            .remove(0);
        assert_eq!(reread.sequence.alphabet(), Alphabet::Protein);
        assert_eq!(reread.sequence.as_str(), "MKL");
        let mut oversized = protein;
        oversized.name = "x".repeat(1025);
        assert!(prepare(&oversized, Format::Fasta).is_err());
    }
}
