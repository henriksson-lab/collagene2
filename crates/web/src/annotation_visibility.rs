//! Display preferences bound to exact retained feature indices, never sequence edits.
use collagene_core::{
    run_lineage::EventBinding,
    workspace::{Workspace, WorkspaceAttachment},
};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

pub const KIND: &str = "annotation-display-v1";
const MAX_BYTES: usize = 2 * 1024 * 1024;
const MAX_HIDDEN: usize = collagene_view_layout::MAX_ITEMS;
type Assets = BTreeMap<String, WorkspaceAttachment>;

#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Record {
    pub schema_version: u16,
    pub source: EventBinding,
    pub feature_count: usize,
    #[serde(deserialize_with = "hidden_indices")]
    pub hidden: BTreeSet<usize>,
}

fn hidden_indices<'de, D: serde::Deserializer<'de>>(d: D) -> Result<BTreeSet<usize>, D::Error> {
    struct Indices;
    impl<'de> serde::de::Visitor<'de> for Indices {
        type Value = BTreeSet<usize>;
        fn expecting(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
            f.write_str("bounded unique annotation indices")
        }
        fn visit_seq<A: serde::de::SeqAccess<'de>>(
            self,
            mut seq: A,
        ) -> Result<Self::Value, A::Error> {
            let mut indices = BTreeSet::new();
            for _ in 0..MAX_HIDDEN {
                let Some(index) = seq.next_element()? else {
                    return Ok(indices);
                };
                if !indices.insert(index) {
                    return Err(serde::de::Error::custom("duplicate hidden index"));
                }
            }
            if seq.next_element::<serde::de::IgnoredAny>()?.is_some() {
                return Err(serde::de::Error::custom("too many hidden annotations"));
            }
            Ok(indices)
        }
    }
    d.deserialize_seq(Indices)
}

impl Record {
    fn validate(&self) -> Result<(), String> {
        if self.schema_version != 1
            || self.hidden.len() > MAX_HIDDEN
            || self.hidden.iter().any(|i| *i >= self.feature_count)
            || [&self.source.document_id, &self.source.event_id]
                .iter()
                .any(|s| s.is_empty() || s.len() > 256 || s.chars().any(char::is_control))
            || self.source.revision.len() != 64
            || !self
                .source
                .revision
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        {
            return Err("Invalid annotation display preferences.".into());
        }
        Ok(())
    }
}

fn records(assets: &Assets) -> Result<Vec<(String, Record)>, String> {
    if assets.len() > collagene_core::workspace::MAX_ATTACHMENTS {
        return Err("Too many workspace attachments.".into());
    }
    let mut result = Vec::new();
    let mut identities = BTreeSet::new();
    for (hash, asset) in assets.iter().filter(|(_, a)| a.kind == KIND) {
        if asset.media_type != "application/json"
            || asset.name != "annotation-display.json"
            || asset.bytes.len() > MAX_BYTES
            || blake3::hash(&asset.bytes).to_hex().as_str() != hash
        {
            return Err("Annotation display attachment identity or metadata differs.".into());
        }
        let record: Record = serde_json::from_slice(&asset.bytes).map_err(|e| e.to_string())?;
        record.validate()?;
        let identity = (
            record.source.document_id.clone(),
            record.source.event_id.clone(),
        );
        if !identities.insert(identity) {
            return Err("Duplicate annotation display preferences for one event.".into());
        }
        result.push((hash.clone(), record));
    }
    Ok(result)
}

pub fn validate_workspace(frame: &Workspace) -> Result<(), String> {
    for (_, record) in records(&frame.attachments)? {
        let event = frame
            .history
            .iter()
            .find(|event| event.id == record.source.event_id)
            .ok_or("Annotation display source event is unavailable.")?;
        let molecule = frame
            .revisions
            .get(&record.source.revision)
            .ok_or("Annotation display source revision is unavailable.")?;
        if event.document_id != record.source.document_id
            || event.revision != record.source.revision
            || molecule.features.len() != record.feature_count
        {
            return Err(
                "Annotation display preferences differ from their exact source event.".into(),
            );
        }
    }
    Ok(())
}

pub fn load(
    assets: &Assets,
    source: &EventBinding,
    count: usize,
) -> Result<BTreeSet<usize>, String> {
    for (_, record) in records(assets)? {
        if record.source == *source {
            if record.feature_count != count {
                return Err("Annotation display feature count changed.".into());
            }
            return Ok(record.hidden);
        }
    }
    Ok(BTreeSet::new())
}

/// Stage only attachments. The host validates the resulting workspace before installation.
pub fn store(assets: &Assets, record: Record) -> Result<Assets, String> {
    record.validate()?;
    let existing = records(assets)?;
    let mut staged = assets.clone();
    for (hash, old) in existing {
        if old.source == record.source {
            staged.remove(&hash);
        }
    }
    if !record.hidden.is_empty() {
        let bytes = serde_json::to_vec(&record).map_err(|e| e.to_string())?;
        if bytes.len() > MAX_BYTES {
            return Err("Annotation display preferences exceed their byte limit.".into());
        }
        let hash = blake3::hash(&bytes).to_hex().to_string();
        if staged.contains_key(&hash) {
            return Err("Annotation display attachment collides with another asset.".into());
        }
        if staged.len() >= collagene_core::workspace::MAX_ATTACHMENTS {
            return Err("Workspace attachment limit reached.".into());
        }
        staged.insert(
            hash,
            WorkspaceAttachment {
                kind: KIND.into(),
                name: "annotation-display.json".into(),
                media_type: "application/json".into(),
                bytes,
            },
        );
    }
    Ok(staged)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn record() -> Record {
        Record {
            schema_version: 1,
            source: EventBinding {
                document_id: "doc".into(),
                event_id: "event".into(),
                revision: "a".repeat(64),
            },
            feature_count: 3,
            hidden: BTreeSet::from([0, 2]),
        }
    }
    #[test]
    fn exact_event_preferences_replace_and_clear_without_aliasing_other_events() {
        let original = record();
        let assets = store(&Assets::new(), original.clone()).unwrap();
        assert_eq!(load(&assets, &original.source, 3).unwrap(), original.hidden);
        let mut later = original.clone();
        later.source.event_id = "later".into();
        assert!(load(&assets, &later.source, 3).unwrap().is_empty());
        let mut changed = original.clone();
        changed.hidden = BTreeSet::from([1]);
        let updated = store(&assets, changed.clone()).unwrap();
        assert_eq!(updated.len(), 1);
        assert_eq!(load(&updated, &original.source, 3).unwrap(), changed.hidden);
        changed.hidden.clear();
        assert!(store(&updated, changed).unwrap().is_empty());
        assert_eq!(assets.len(), 1);
    }
    #[test]
    fn malformed_or_tampered_preferences_cannot_change_the_input_assets() {
        let original = record();
        let assets = store(&Assets::new(), original.clone()).unwrap();
        assert!(load(&assets, &original.source, 4).is_err());
        let mut invalid = original.clone();
        invalid.hidden.insert(3);
        assert!(store(&assets, invalid).is_err());
        let mut tampered = assets.clone();
        tampered.values_mut().next().unwrap().bytes.push(b' ');
        assert!(load(&tampered, &original.source, 3).is_err());
        let mut json = serde_json::to_value(&original).unwrap();
        json["hidden"] = serde_json::json!([0, 0]);
        assert!(serde_json::from_value::<Record>(json).is_err());
        assert_eq!(load(&assets, &original.source, 3).unwrap(), original.hidden);
    }
    #[test]
    fn native_archive_preferences_require_the_exact_owner_event_and_feature_count() {
        use collagene_model::{Feature, Location, Molecule, Topology};
        let mut molecule = Molecule::dna("source", "ACGT", Topology::Linear).unwrap();
        molecule.features.push(Feature {
            kind: "gene".into(),
            location: Location::range(0, 3),
            qualifiers: vec![],
            needs_review: false,
        });
        let mut frame =
            Workspace::from_documents("display", vec![molecule.clone(), molecule], 0).unwrap();
        let source = EventBinding {
            document_id: frame.documents[0].id.clone(),
            event_id: frame.history[0].id.clone(),
            revision: frame.documents[0].head.clone(),
        };
        let preference = Record {
            schema_version: 1,
            source,
            feature_count: 1,
            hidden: BTreeSet::from([0]),
        };
        frame.attachments = store(&frame.attachments, preference.clone()).unwrap();
        validate_workspace(&frame).unwrap();
        let restored = collagene_core::workspace::load_workspace(
            &collagene_core::workspace::save_workspace(&frame).unwrap(),
        )
        .unwrap();
        validate_workspace(&restored).unwrap();
        assert_eq!(
            load(&restored.attachments, &preference.source, 1).unwrap(),
            preference.hidden
        );
        let mut alias = preference.clone();
        alias.source.document_id = frame.documents[1].id.clone();
        frame.attachments = store(&Assets::new(), alias).unwrap();
        assert!(validate_workspace(&frame).is_err());
        let mut wrong_count = preference;
        wrong_count.feature_count = 2;
        frame.attachments = store(&Assets::new(), wrong_count).unwrap();
        assert!(validate_workspace(&frame).is_err());
    }
}
