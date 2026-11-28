use collagene_core::{
    documents::workspace::{self, Workspace, WorkspaceAttachment},
    protein::structure::{self as structure_view, coordinates::Format},
};
use collagene_library_sharing::{self as sharing, Selection};
use collagene_model::{Molecule, Topology};
const PDB:&[u8]=b"REMARK   1 UNLINKED PRIVATE COORDINATE SOURCE\r\nATOM      1  CA  ALA A   1      11.104  13.207   2.100  1.00 20.00           C  \r\nEND\r\n";
const CIF:&[u8]=b"data_unlinked\n# UNLINKED PRIVATE COORDINATE SOURCE\nloop_\n_atom_site.group_PDB\n_atom_site.id\n_atom_site.type_symbol\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\nATOM 1 C CA ALA A 1 11.104 13.207 2.100\n";
fn fixture() -> (Workspace, Selection) {
    let w = Workspace::from_documents(
        "Private full workspace",
        vec![
            Molecule::dna("Shared plasmid", "ACGT", Topology::Circular).unwrap(),
            Molecule::dna("PRIVATE other plasmid", "GGAA", Topology::Circular).unwrap(),
        ],
        0,
    )
    .unwrap();
    let d = &w.documents[0];
    let selection = Selection {
        document_id: d.id.clone(),
        revision: d.head.clone(),
        head_event: d.head_event.clone(),
    };
    (w, selection)
}
#[test]
fn unlinked_originals_survive_full_backup_but_do_not_block_or_expand_selected_publication() {
    let (mut w, selection) = fixture();
    let mut hashes = vec![];
    for (bytes, format, name) in [
        (PDB, Format::Pdb, "unlinked.pdb"),
        (CIF, Format::Mmcif, "unlinked.cif"),
    ] {
        // Authored sources really are importable; publication does not need to
        // reinterpret their atom data or guess sequence correspondence.
        let prepared = structure_view::prepare(bytes, format, None).unwrap();
        assert_eq!(prepared.atoms.len(), 1);
        hashes.push(
            structure_view::retain(&mut w.attachments, name, bytes.to_vec(), format).unwrap(),
        );
    }
    let before = serde_json::to_value(&w).unwrap();
    let loaded = workspace::load_workspace(&workspace::save_workspace(&w).unwrap()).unwrap();
    assert_eq!(loaded.attachments[&hashes[0]].bytes, PDB);
    assert_eq!(loaded.attachments[&hashes[1]].bytes, CIF);
    let inventory = sharing::inspect_document_closure(&w, &selection).unwrap();
    assert_eq!(inventory.engine_version, "library-document-closure-4");
    assert!(inventory.attachments.is_empty());
    assert_eq!(inventory.documents.len(), 1);
    assert_eq!(inventory.omitted.attachments, 2);
    let published = sharing::document_closure(&w, &selection).unwrap();
    let checked = workspace::load_workspace(&published.archive).unwrap();
    assert_eq!(checked.documents.len(), 1);
    assert!(checked.attachments.is_empty());
    assert_eq!(
        checked.head(&selection.document_id).unwrap().name,
        "Shared plasmid"
    );
    let private_view = sharing::project_private_document_closure(&w, &selection).unwrap();
    assert!(private_view.workspace.attachments.is_empty());
    assert_eq!(private_view.workspace.documents.len(), 1);
    assert_eq!(serde_json::to_value(&w).unwrap(), before);
}
#[test]
fn explicit_coordinate_omission_does_not_turn_unknown_attachment_schemas_into_workspace_only_data()
{
    let (mut w, selection) = fixture();
    for kind in ["structure-source-pdb-v0", "future-foreign-key-contract"] {
        let bytes = b"unknown ownership contract".to_vec();
        let hash = blake3::hash(&bytes).to_hex().to_string();
        w.attachments.insert(
            hash,
            WorkspaceAttachment {
                kind: kind.into(),
                name: "Unknown reference data".into(),
                media_type: "application/octet-stream".into(),
                bytes,
            },
        );
        assert!(
            matches!(sharing::inspect_document_closure(&w,&selection),Err(sharing::Error::UnsupportedAttachment(k)) if k==kind)
        );
        w.attachments.clear();
    }
}
