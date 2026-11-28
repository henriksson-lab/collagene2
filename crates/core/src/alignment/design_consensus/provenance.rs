//! Bounded retained-original facts, not authentication or sequencing quality.
use super::*;
use base64::{Engine, engine::general_purpose::STANDARD};
use collagene_model::Molecule;
use std::borrow::Cow;
const MAX_ORIGINAL: usize = 4 * 1024 * 1024;
fn inspect(molecule: &Molecule) -> OriginalEvidence {
    let raw = if let Some(raw) = &molecule.original_genbank {
        Some(("GenBank", Cow::Borrowed(raw.as_bytes()), None))
    } else if let Some(raw) = molecule.metadata.get("embl.original_record") {
        Some(("EMBL", Cow::Borrowed(raw.as_bytes()), None))
    } else if let Some(encoded) = molecule.metadata.get("source:original_fasta_base64") {
        if encoded.len() > MAX_ORIGINAL.div_ceil(3) * 4 {
            return invalid(
                "FASTA",
                "Retained original encoding exceeds the bounded evidence limit",
            );
        }
        match STANDARD.decode(encoded) {
            Ok(bytes) => Some((
                "FASTA",
                Cow::Owned(bytes),
                molecule.metadata.get("source:original_fasta_blake3"),
            )),
            Err(_) => return invalid("FASTA", "Retained original FASTA encoding is invalid"),
        }
    } else {
        None
    };
    let Some((format, bytes, declared)) = raw else {
        return OriginalEvidence{state:OriginalState::Unavailable,format:None,blake3:None,byte_count:0,
        message:"No supported original record bytes were retained. Exact normalized event calls remain available; provider identity and read quality are not inferred.".into()};
    };
    if bytes.len() > MAX_ORIGINAL {
        return invalid(
            format,
            "Retained original record exceeds the bounded evidence limit",
        );
    }
    let hash = blake3::hash(&bytes).to_hex().to_string();
    if declared.is_some_and(|declared| declared != &hash) {
        return invalid(
            format,
            "Retained original record digest differs from its declared hash",
        );
    }
    let parsed = crate::formats::sequence::import_bytes(&bytes);
    let state = match parsed {
        Ok(parsed) if parsed.molecules.len() == 1 => {
            if parsed.molecules[0].sequence == molecule.sequence {
                OriginalState::IntegrityCheckedMatchesCalls
            } else {
                OriginalState::IntegrityCheckedDifferentCalls
            }
        }
        _ => OriginalState::InvalidOrUnsupported,
    };
    let message = match state {
        OriginalState::IntegrityCheckedMatchesCalls => {
            "Retained bytes reparsed to the selected current nucleotide calls. This does not authenticate provider/sample identity or establish experimental confidence."
        }
        OriginalState::IntegrityCheckedDifferentCalls => {
            "Retained original bytes and current calls differ; historical file evidence was not substituted for the selected immutable event."
        }
        _ => {
            "Retained original bytes could not be checked as one supported sequence record; no original/current correspondence claim is made."
        }
    };
    OriginalEvidence {
        state,
        format: Some(format.into()),
        blake3: Some(hash),
        byte_count: bytes.len(),
        message: message.into(),
    }
}
fn invalid(format: &str, message: &str) -> OriginalEvidence {
    OriginalEvidence {
        state: OriginalState::InvalidOrUnsupported,
        format: Some(format.into()),
        blake3: None,
        byte_count: 0,
        message: message.into(),
    }
}
pub(super) fn capture(
    role: SourceRole,
    event: &crate::material::lineage::EventBinding,
    molecule: &Molecule,
) -> SourceProvenance {
    // Known documentary fields only. URL/path/credential strings are not copied
    // into receipts, and this owner never contacts a provider.
    let fields = [
        "source:provider",
        "source:format",
        "source:identifier",
        "source:accession_version",
        "source:retrieved_at",
        "source:record_index",
        "source:content_blake3",
    ]
    .iter()
    .filter_map(|key| {
        molecule.metadata.get(*key).map(|value| {
            let mut bytes = 0usize;
            let displayed: String = value
                .chars()
                .filter(|c| !c.is_control())
                .take_while(|c| {
                    bytes += c.len_utf8();
                    bytes <= 2048
                })
                .collect();
            ProvenanceField {
                key: (*key).into(),
                value_complete: displayed == *value,
                value: displayed,
            }
        })
    })
    .filter(|field| !field.value.trim().is_empty())
    .collect();
    SourceProvenance {
        role,
        event: event.clone(),
        topology: molecule.topology,
        fields,
        original: inspect(molecule),
    }
}
