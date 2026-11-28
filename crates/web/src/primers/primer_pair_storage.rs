//! Unfinished editable text, never a saved binding or a positive PCR result.
use crate::primers::primer_pair_model::PairDraft;
use crate::primers::primer_pair_model::anchor_reference;
use crate::primers::primer_pair_model::selected_guard;
#[cfg(test)]
use collagene_core::documents::reference_edit::Operation as ReferenceEdit;
use collagene_core::material::lineage::EventBinding;
use collagene_model::Molecule;
use serde::{Deserialize, Serialize};

pub const ATTACHMENT_KIND: &str = "primer-pair-draft-v2";
const MAX_BYTES: usize = 64 * 1024;

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Saved {
    version: u16,
    draft: PairDraft,
}

fn text(value: &str, chars: usize) -> Result<(), String> {
    if value.len() > chars * 4 || value.chars().count() > chars {
        return Err("Primer pair draft exceeds its editable field limits.".into());
    }
    Ok(())
}
fn guard(draft: &PairDraft) -> Result<(), String> {
    for oligo in [&draft.forward, &draft.reverse] {
        // Incomplete or currently invalid names/DNA are deliberate editable text.
        // Scientific conversion still uses OligoDraft::primer after restoration.
        text(&oligo.name, 256)?;
        text(&oligo.annealing, 1600)?;
        text(&oligo.tail, 1600)?;
        if let Some(a) = &oligo.anchor {
            selected_guard(&a.source, a.start, a.end, &a.reference)?;
        }
    }
    Ok(())
}

pub fn encode(draft: &PairDraft) -> Result<Vec<u8>, String> {
    guard(draft)?;
    let bytes = serde_json::to_vec(&Saved {
        version: 2,
        draft: draft.clone(),
    })
    .map_err(|e| e.to_string())?;
    if bytes.len() > MAX_BYTES {
        return Err("Primer pair draft exceeds64KiB.".into());
    }
    Ok(bytes)
}

/// The host resolves exact retained events, including historical source heads.
/// Equal DNA in a different document/event cannot stand in for that event.
pub fn decode(
    bytes: &[u8],
    mut resolve: impl FnMut(&EventBinding) -> Result<Molecule, String>,
) -> Result<PairDraft, String> {
    if bytes.len() > MAX_BYTES {
        return Err("Primer pair draft exceeds64KiB.".into());
    }
    let saved: Saved = serde_json::from_slice(bytes).map_err(|e| e.to_string())?;
    if saved.version != 2 {
        return Err("Unsupported primer pair draft version.".into());
    }
    guard(&saved.draft)?;
    for oligo in [&saved.draft.forward, &saved.draft.reverse] {
        let Some(a) = &oligo.anchor else { continue };
        let source = resolve(&a.source)?;

        let expected = anchor_reference(&source, &a.source, a.start, a.end, a.reverse)?;
        if a.reference != expected {
            return Err(
                "Primer selection anchor bases differ from the actual retained event.".into(),
            );
        }
    }
    Ok(saved.draft)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primers::primer_pair_model::SelectionSeed;
    use collagene_core::documents::history::DocumentHistory;
    use collagene_model::Topology;

    fn binding(doc: &DocumentHistory) -> EventBinding {
        EventBinding {
            document_id: doc.id.clone(),
            event_id: doc.head_event().into(),
            revision: doc.head_revision().into(),
        }
    }

    #[test]
    fn unfinished_invalid_text_round_trips_without_becoming_a_primer() {
        let mut draft = PairDraft::default();
        draft.forward.name = " α primer ".into();
        draft.forward.annealing = "ac g\n?".into();
        draft.forward.tail = "🧬".into();
        draft.reverse.name.clear();
        let bytes = encode(&draft).unwrap();
        let decoded = decode(&bytes, |_| panic!("unanchored text needs no source")).unwrap();
        assert_eq!(decoded, draft);
        assert!(decoded.forward.primer().is_err());
        assert!(decoded.reverse.primer().is_err());
    }

    #[test]
    fn original_historical_event_and_reverse_selected_bases_survive_archive_exchange() {
        let molecule = Molecule::dna("original", "AAGTCGATCC", Topology::Linear).unwrap();
        let mut doc = DocumentHistory::new("template", molecule.clone()).unwrap();
        let exact = binding(&doc);
        let mut draft = PairDraft::default();
        draft.reverse = SelectionSeed::capture(&molecule, &exact, 1, 6)
            .unwrap()
            .oligo("retained reverse", true, "ATG")
            .unwrap();
        assert_eq!(draft.reverse.anchor.as_ref().unwrap().reference, "CGACT");
        // An intentional mismatch stays in the draft, independently of its anchor.
        draft.reverse.annealing = "CGGCT".into();
        doc.commit(
            vec![
                ReferenceEdit::Rename {
                    name: "edited head".into(),
                },
                ReferenceEdit::ReplaceSequence {
                    start: 0,
                    end: 10,
                    replacement: "TTTTTTTTTT".into(),
                },
            ],
            "edit",
        )
        .unwrap();
        let mut workspace = collagene_core::documents::workspace::Workspace::from_documents(
            "pair",
            vec![doc.molecule.clone()],
            0,
        )
        .unwrap();
        doc.contribute(&mut workspace).unwrap();
        workspace.active_document = doc.id.clone();
        let bytes = encode(&draft).unwrap();
        workspace.attachments.insert(
            blake3::hash(&bytes).to_hex().to_string(),
            collagene_core::documents::workspace::WorkspaceAttachment {
                kind: ATTACHMENT_KIND.into(),
                name: "pair.json".into(),
                media_type: "application/json".into(),
                bytes: bytes.clone(),
            },
        );
        let archive = collagene_core::documents::workspace::save_workspace(&workspace).unwrap();
        let retained = collagene_core::documents::workspace::load_workspace(&archive).unwrap();
        let sources = DocumentHistory::restore_all(&retained).unwrap();
        let saved_asset = retained
            .attachments
            .values()
            .find(|a| a.kind == ATTACHMENT_KIND)
            .unwrap();
        assert_eq!(saved_asset.bytes, bytes);
        let decoded = decode(&saved_asset.bytes, |b| {
            sources
                .iter()
                .find(|d| d.id == b.document_id)
                .ok_or("missing exact document")?
                .revision_molecule(&b.event_id)
                .map_err(|e| e.to_string())
        })
        .unwrap();
        assert_eq!(decoded, draft);
        assert_eq!(
            decoded.reverse.comparison().unwrap().unwrap().mismatches,
            vec![2]
        );
        // Same sequence/revision in an unrelated document cannot supply the event.
        let alias = DocumentHistory::new("alias", molecule).unwrap();
        assert!(
            decode(&bytes, |b| {
                if alias.id != b.document_id {
                    return Err("wrong document".into());
                }
                alias
                    .revision_molecule(&b.event_id)
                    .map_err(|e| e.to_string())
            })
            .is_err()
        );
    }

    #[test]
    fn unsupported_malformed_stale_and_oversized_drafts_fail_before_restore() {
        let molecule = Molecule::dna("source", "AAGTCGATCC", Topology::Linear).unwrap();
        let doc = DocumentHistory::new("template", molecule.clone()).unwrap();
        let mut draft = PairDraft::default();
        draft.forward = SelectionSeed::capture(&molecule, &binding(&doc), 0, 5)
            .unwrap()
            .oligo("f", false, "")
            .unwrap();
        let bytes = encode(&draft).unwrap();
        let resolve = |_: &EventBinding| Ok(molecule.clone());
        let mut value: serde_json::Value = serde_json::from_slice(&bytes).unwrap();
        value["version"] = 3.into();
        assert!(decode(&serde_json::to_vec(&value).unwrap(), resolve).is_err());
        value["version"] = 2.into();
        value["draft"]["forward"]["anchor"]["reference"] = "TTTTT".into();
        assert!(decode(&serde_json::to_vec(&value).unwrap(), resolve).is_err());
        assert!(
            decode(&bytes, |_| Ok(Molecule::dna(
                "wrong",
                "TTTTTTTTTT",
                Topology::Linear
            )
            .unwrap()))
            .is_err()
        );
        assert!(
            decode(&vec![b' '; MAX_BYTES + 1], |_| panic!(
                "oversize must fail first"
            ))
            .is_err()
        );
        draft.reverse.annealing = "A".repeat(1601);
        assert!(encode(&draft).is_err());
    }

    #[test]
    fn wrapped_anchor_reopens_exact_historical_event_and_rejects_unsupported_schema() {
        let source = Molecule::dna("original circle", "ACGTAAGC", Topology::Circular).unwrap();
        let mut doc = DocumentHistory::new("circle", source.clone()).unwrap();
        let exact = binding(&doc);
        let mut draft = PairDraft::default();
        draft.reverse = SelectionSeed::capture(&source, &exact, 6, 11)
            .unwrap()
            .oligo("wrapped reverse", true, "GG")
            .unwrap();
        draft.reverse.annealing = "incomplete ?".into();
        let bytes = encode(&draft).unwrap();
        doc.commit(
            vec![
                ReferenceEdit::Rename {
                    name: "later circle".into(),
                },
                ReferenceEdit::ReplaceSequence {
                    start: 0,
                    end: 8,
                    replacement: "TTTTTTTT".into(),
                },
            ],
            "head changed",
        )
        .unwrap();
        let restored = decode(&bytes, |binding| {
            doc.revision_molecule(&binding.event_id)
                .map_err(|e| e.to_string())
        })
        .unwrap();
        assert_eq!(restored, draft);
        assert_eq!(restored.reverse.anchor.as_ref().unwrap().reference, "CGTGC");
        assert!(restored.reverse.primer().is_err());
        assert!(decode(&bytes, |_| Ok(doc.molecule.clone())).is_err());
        let mut legacy: serde_json::Value = serde_json::from_slice(&bytes).unwrap();
        legacy["version"] = 1.into();
        assert!(
            decode(
                &serde_json::to_vec(&legacy).unwrap(),
                |_| Ok(source.clone())
            )
            .is_err()
        );
    }
    #[test]
    fn current_pair_requires_explicit_anchor_and_rejects_unknown_fields() {
        let bytes = encode(&PairDraft::default()).unwrap();
        let current: serde_json::Value = serde_json::from_slice(&bytes).unwrap();
        let mut missing = current.clone();
        missing["draft"]["forward"]
            .as_object_mut()
            .unwrap()
            .remove("anchor");
        assert!(
            decode(&serde_json::to_vec(&missing).unwrap(), |_| panic!(
                "No anchor may resolve"
            ))
            .is_err()
        );
        let mut unknown = current;
        unknown["draft"]["obsolete"] = true.into();
        assert!(
            decode(&serde_json::to_vec(&unknown).unwrap(), |_| panic!(
                "No anchor may resolve"
            ))
            .is_err()
        );
    }
}
