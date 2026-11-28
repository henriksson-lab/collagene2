//! Declared CDS coordinate/code pins, using canonical extraction and translation.
use super::*;
use crate::{construct_qc::AnnotationAssumption, sequence::translation};
use collagene_model::{
    Alphabet, Location, Molecule,
    annotations::biological_path::{BiologicalCall, BiologicalPath, PathLimits},
};

pub(super) struct Coding {
    pub path: BiologicalPath,
    pub offset: usize,
    pub source_length: usize,
    pub topology: collagene_model::Topology,
    pub oriented: String,
    pub translated: translation::CdsTranslation,
}
impl Coding {
    pub fn load(
        molecule: &Molecule,
        event_revision: &str,
        assumption: &AnnotationAssumption,
        code: u8,
        remaining: &mut usize,
    ) -> Result<Self, Error> {
        if molecule.sequence.alphabet() == Alphabet::Protein {
            return Err(Error::Unsupported(
                "Native CDS comparison requires nucleotide data".into(),
            ));
        }
        if assumption.source_revision != event_revision {
            return Err(Error::Unsupported(
                "CDS declaration belongs to another source revision".into(),
            ));
        }
        let feature = molecule
            .features
            .get(assumption.feature_index)
            .ok_or_else(|| Error::Unsupported("Declared CDS annotation is unavailable".into()))?;
        collagene_model::bounds::feature_preflight(&[std::slice::from_ref(feature)])?;
        crate::encoding::bounded_json::preflight(feature, 256 << 10)?;
        if feature.needs_review || !feature.kind.eq_ignore_ascii_case("CDS") {
            return Err(Error::Unsupported("CDS role is missing or still requires annotation review; a declaration cannot silently verify it".into()));
        }
        if crate::construct_qc::annotation_revision(feature)? != assumption.annotation_revision {
            return Err(Error::Unsupported(
                "CDS declaration does not pin the actual annotation".into(),
            ));
        }
        let path = BiologicalPath::from_location(
            &feature.location,
            molecule.sequence.len(),
            PathLimits {
                max_calls: 2_000_000,
                max_spans: 4096,
            },
        )?;
        if path.len() > *remaining {
            return Err(Error::TooLarge);
        }
        *remaining -= path.len();
        let mut effective = feature.clone();
        if !effective
            .qualifiers
            .iter()
            .any(|(key, _)| key == "transl_table")
        {
            effective
                .qualifiers
                .push(("transl_table".into(), Some(code.to_string())));
        }
        let translated =
            translation::translate_cds_validated(molecule, &effective, event_revision)?;
        if translated.translation.options.genetic_code != code {
            return Err(Error::Unsupported(
                "Selected genetic code conflicts with the annotation's exact table".into(),
            ));
        }
        let oriented = feature.location.extract(&molecule.sequence)?;
        if oriented.len() != path.len() {
            return Err(Error::Unsupported(
                "Canonical biological path and native CDS extraction disagree".into(),
            ));
        }
        Ok(Self {
            path,
            offset: usize::from(translated.codon_start - 1),
            source_length: molecule.sequence.len(),
            topology: molecule.topology,
            oriented,
            translated,
        })
    }
    pub fn coding_len(&self) -> usize {
        self.path.len().saturating_sub(self.offset)
    }
    pub fn codon(&self, index: usize) -> Option<&str> {
        let start = index.checked_mul(3)?.checked_add(self.offset)?;
        self.oriented.get(start..start.checked_add(3)?)
    }
    pub fn codon_coordinates(&self, index: usize) -> Option<Location> {
        let start = index.checked_mul(3)?.checked_add(self.offset)?;
        let calls = (0..3)
            .map(|offset| self.path.call_at(start + offset))
            .collect::<Option<Vec<_>>>()?;
        Some(calls_location(&calls))
    }
    pub fn amino(&self, index: usize) -> Option<char> {
        self.translated
            .translation
            .protein
            .as_bytes()
            .get(index)
            .copied()
            .map(char::from)
    }
}
pub(super) fn calls_location(calls: &[BiologicalCall]) -> Location {
    let parts = calls
        .iter()
        .map(|call| {
            let location = Location::range(call.position, call.position + 1);
            if call.reverse {
                Location::Complement(Box::new(location))
            } else {
                location
            }
        })
        .collect::<Vec<_>>();
    if parts.len() == 1 {
        parts.into_iter().next().expect("one biological call")
    } else {
        Location::Join(parts)
    }
}
