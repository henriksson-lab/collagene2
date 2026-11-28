//! WebMCP access to the live declaration, without serialized source bodies or receipts.
use crate::assembly_draft_model::Draft;
use serde::{Deserialize, Serialize};
use serde_json::Value;
use std::io::Write;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct Pin {
    pub workspace_epoch: u64,
    pub draft_revision: String,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct ReviewToken {
    pub owner: u64,
    pub generation: u32,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct ReviewRequest {
    pub expected: Pin,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct CancelRequest {
    pub expected: Pin,
    pub review: ReviewToken,
}

struct Hash {
    inner: blake3::Hasher,
    bytes: usize,
}
impl Write for Hash {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(bytes.len())
            .filter(|n| *n <= 4 * 1024 * 1024)
            .ok_or_else(|| std::io::Error::other("Assembly tool declaration exceeds 4 MiB"))?;
        self.inner.update(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn slots(draft: &Draft) -> Vec<Value> {
    draft
        .slots
        .iter()
        .map(|slot| {
            serde_json::json!({
                "label":slot.label, "source":slot.captured.as_ref().map(|s|&s.source),
                "selector":slot.selector, "digest_geometry":slot.digest_geometry,
                "material_selection":slot.material_origin.as_ref().map(|origin|origin.selection()),
                "material_revision":slot.material_origin.as_ref().map(|origin|origin.revision()),
                "bottom_source":slot.anneal_bottom.as_ref().map(|s|&s.source)
            })
        })
        .collect()
}
fn methods(draft: &Draft) -> crate::cloning_methods_model::MethodsState {
    let mut methods = draft.methods.clone();
    for row in &mut methods.rows {
        row.run_origin = None;
    }
    methods
}
fn row_origins(
    draft: &Draft,
) -> Vec<Option<(&str, collagene_core::run_material::MaterialSelection)>> {
    draft
        .methods
        .rows
        .iter()
        .map(|row| {
            row.run_origin
                .as_ref()
                .map(|origin| (origin.revision(), origin.selection()))
        })
        .collect()
}
/// Declaration key excludes immutable source DNA/capsules; their captured exact
/// source revisions and cached origin identities bind them without replay/copy.
pub(crate) fn pin(draft: &Draft, epoch: u64) -> Result<Pin, String> {
    if draft.slots.len() > 12 {
        return Err("Assembly tools support at most 12 declared slots.".into());
    }
    let alternatives = draft
        .alternatives
        .as_ref()
        .map(|_| crate::assembly_batch_projection::declaration_revision(draft))
        .transpose()?;
    let mut hash = Hash {
        inner: blake3::Hasher::new(),
        bytes: 0,
    };
    serde_json::to_writer(
        &mut hash,
        &(
            "collagene-webmcp-assembly",
            collagene_core::verification::ENGINE_VERSION,
            draft.version,
            &draft.name,
            draft.circular,
            slots(draft),
            methods(draft),
            row_origins(draft),
            &draft.junction_edits,
            &draft.substitutions,
            &draft.edits,
            alternatives,
        ),
    )
    .map_err(|e| e.to_string())?;
    Ok(Pin {
        workspace_epoch: epoch,
        draft_revision: hash.inner.finalize().to_hex().to_string(),
    })
}
pub(crate) fn require_pin(draft: &Draft, epoch: u64, expected: &Pin) -> Result<(), String> {
    if pin(draft, epoch)? != *expected {
        return Err(
            "The assembly draft or workspace changed. Read its current context again.".into(),
        );
    }
    Ok(())
}
pub(crate) fn context(draft: &Draft, epoch: u64) -> Result<Value, String> {
    Ok(serde_json::json!({
        "scope":"active_assembly_draft", "pin":pin(draft,epoch)?,
        "name":draft.name,"circular":draft.circular,"active_slot":draft.active_slot,
        "slots":slots(draft),"method_settings_without_material_capsules":methods(draft),"preparation_material_origins":row_origins(draft),
        "junction_edits":draft.junction_edits,"substitutions":draft.substitutions,"edits":draft.edits,
        "active_library_construct":crate::assembly_alternatives::identity(draft),
        "notice":"Pinned declared requirements and unfinished settings. Intended-reference review is separate from method verification. Primer proposals and application remain in the assembly review controls."
    }))
}

#[cfg(target_arch = "wasm32")]
pub(crate) mod bridge {
    use super::*;
    use leptos::prelude::*;
    use std::sync::{Arc, atomic::AtomicBool};
    #[derive(Clone)]
    pub(crate) enum Command {
        Status(Pin),
        Review(Pin),
        Cancel(Pin, ReviewToken),
        Propose(crate::web_mcp_design::Request),
        Discard(crate::web_mcp_design::DiscardRequest),
    }
    #[derive(Clone)]
    pub(crate) struct Endpoint {
        pub alive: Arc<AtomicBool>,
        pub command: Callback<Command, Result<Value, String>>,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn source_event_selector_and_physical_settings_are_exact_but_bodies_are_not_exposed() {
        use crate::{assembly_draft_model::CapturedSource, source_identity::SourceInfo};
        use collagene_core::history::DocumentHistory;
        use collagene_model::{Molecule, Topology};
        let mut draft = Draft::default();
        let doc = DocumentHistory::new(
            "original",
            Molecule::dna("Original source", "ACGTACGT", Topology::Linear).unwrap(),
        )
        .unwrap();
        let source = SourceInfo {
            binding: collagene_core::run_lineage::EventBinding {
                document_id: doc.id.clone(),
                event_id: doc.head_event().into(),
                revision: doc.head_revision().into(),
            },
            name: doc.molecule.name.clone(),
            length: 8,
            alphabet: doc.molecule.sequence.alphabet(),
            topology: doc.molecule.topology,
        };
        draft.slots[0].captured = Some(std::sync::Arc::new(CapturedSource {
            source,
            molecule: doc.molecule,
        }));
        draft.slots[0].selector = collagene_core::assembly_intent::Selector::ReferenceRange {
            start: 0,
            end: 8,
            reverse: false,
        };
        let expected = pin(&draft, 9).unwrap();
        let output = context(&draft, 9).unwrap();
        assert_eq!(
            output["slots"][0]["source"]["binding"]["event_id"],
            draft.slots[0]
                .captured
                .as_ref()
                .unwrap()
                .source
                .binding
                .event_id
        );
        assert!(output["slots"][0].get("molecule").is_none());
        assert!(output["slots"][0].get("sequence").is_none());
        let mut different = draft.clone();
        std::sync::Arc::make_mut(different.slots[0].captured.as_mut().unwrap())
            .source
            .binding
            .event_id
            .push_str("alias");
        assert!(require_pin(&different, 9, &expected).is_err());
        let mut different = draft.clone();
        different.slots[0].selector = collagene_core::assembly_intent::Selector::ReferenceRange {
            start: 1,
            end: 8,
            reverse: false,
        };
        assert!(require_pin(&different, 9, &expected).is_err());
        let mut different = draft.clone();
        different.methods.rows[0].reverse = true;
        assert!(require_pin(&different, 9, &expected).is_err());
        require_pin(&draft, 9, &expected).unwrap();
    }
    #[test]
    fn unfinished_raw_settings_and_workspace_epoch_bind_review_admission() {
        let draft = Draft::default();
        let before = pin(&draft, 4).unwrap();
        let mut changed = draft.clone();
        changed.name.push_str(" edited");
        assert!(require_pin(&changed, 4, &before).is_err());
        let mut changed = draft.clone();
        changed
            .methods
            .homology
            .get_or_insert_with(Default::default)
            .cuts = "unfinished".into();
        assert!(require_pin(&changed, 4, &before).is_err());
        assert!(require_pin(&draft, 5, &before).is_err());
        require_pin(&draft, 4, &before).unwrap();
    }
}
