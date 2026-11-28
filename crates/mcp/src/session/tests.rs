use crate::server::CollageneServer;
use crate::session::Artifact;
use crate::session::ArtifactData;
use crate::session::MAX_ITEMS;
use collagene_core::documents::workspace;
use collagene_model::Molecule;
use serde_json::json;

#[cfg(test)]
mod finalized_workspace_atomic_capacity {
    use super::{Artifact, ArtifactData, CollageneServer, MAX_ITEMS, Molecule, json, workspace};
    #[test]
    fn both_finalization_artifacts_are_reserved_before_any_workspace_publication() {
        let server = CollageneServer::default();
        {
            let mut store = server.session.store.lock().unwrap();
            for index in 0..MAX_ITEMS - 1 {
                store.artifacts.insert(
                    format!("reserved-{index}"),
                    Artifact {
                        data: ArtifactData::Text(String::new()),
                        mime_type: "text/plain",
                        encoding: "utf-8",
                        filename: "reserved.txt",
                    },
                );
            }
        }
        let source = Molecule::dna("source", "ACGT", collagene_model::Topology::Linear).unwrap();
        let workspace =
            workspace::Workspace::from_documents("Atomic fixture", vec![source], 0).unwrap();
        let native = workspace::save_workspace(&workspace).unwrap();
        let native_id = blake3::hash(&native).to_hex().to_string();
        let text = "{\"fixture\":\"receipt\"}".to_string();
        let report_id = blake3::hash(text.as_bytes()).to_hex().to_string();
        assert!(
            server
                .session
                .commit_workspace_outputs(workspace, Some((text, json!({}))), Some(native))
                .is_err()
        );
        let store = server.session.store.lock().unwrap();
        assert_eq!(store.artifacts.len(), MAX_ITEMS - 1);
        assert!(store.workspaces.is_empty());
        assert!(store.documents.is_empty());
        assert_eq!(store.used, 0);
        assert!(!store.artifacts.contains_key(&native_id));
        assert!(!store.artifacts.contains_key(&report_id));
    }
}

#[cfg(test)]
mod structure_artifact_tests {
    use super::{CollageneServer, MAX_ITEMS};
    #[test]
    fn source_and_report_capacity_failure_publishes_neither_artifact() {
        let server = CollageneServer::default();
        for n in 0..MAX_ITEMS - 1 {
            server
                .session
                .text_artifact(
                    format!("Occupied artifact {n}"),
                    "text/plain",
                    "occupied.txt",
                )
                .unwrap();
        }
        let source = "ATOM      1  CA  ALA A   1      11.104  13.207   2.100  1.00 20.00           C  \nEND\n";
        let report = collagene_core::protein::structure::prepare(
            source.as_bytes(),
            collagene_core::protein::structure::coordinates::Format::Pdb,
            None,
        )
        .unwrap();
        let before = server.session.store.lock().unwrap().used;
        assert!(
            server
                .session
                .structure_output(source.into(), &report)
                .is_err()
        );
        let store = server.session.store.lock().unwrap();
        assert_eq!(store.used, before);
        assert_eq!(store.artifacts.len(), MAX_ITEMS - 1);
        assert!(!store.artifacts.contains_key(&report.source_hash));
        assert!(store.documents.is_empty());
        assert!(store.plans.is_empty());
        assert!(store.workspaces.is_empty());
    }
}
