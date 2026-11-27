//! Native credential/permission tests only. No HTTP requests or listener.
use collagene_lab_library::collaboration::{ActorId, Author, Role};
use collagene_library_protocol as wire;
use collagene_server::{AccessKey, MemberCredential, MemberRegistry};
use std::path::PathBuf;
struct Temp(PathBuf);
impl Temp {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "collagene-member-fixture-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let mut b = std::fs::DirBuilder::new();
        #[cfg(unix)]
        {
            use std::os::unix::fs::DirBuilderExt;
            b.mode(0o700);
        }
        b.create(&p).unwrap();
        Self(p)
    }
    fn key(&self, name: &str) -> (AccessKey, String) {
        let p = self.0.join(name);
        AccessKey::generate_file(&p).unwrap();
        let header = format!("Bearer {}", std::fs::read_to_string(&p).unwrap().trim());
        (AccessKey::read(&p).unwrap(), header)
    }
    fn config(&self, content: &str) -> PathBuf {
        use std::io::Write;
        let p = self.0.join("members.json");
        let mut options = std::fs::OpenOptions::new();
        options.write(true).create(true).truncate(true);
        #[cfg(unix)]
        {
            use std::os::unix::fs::OpenOptionsExt;
            options.mode(0o600);
        }
        options
            .open(&p)
            .unwrap()
            .write_all(content.as_bytes())
            .unwrap();
        p
    }
}
impl Drop for Temp {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
fn credential(key: AccessKey, id: &str, role: Role) -> MemberCredential {
    MemberCredential::new(
        key,
        Author::new(ActorId::parse(id).unwrap(), id).unwrap(),
        role,
    )
    .unwrap()
}
#[test]
fn registry_enforces_unique_identities_keys_and_atomic_revocation() {
    let t = Temp::new();
    let (key, header) = t.key("a.key");
    let (other, other_header) = t.key("b.key");
    let registry = MemberRegistry::new(vec![
        credential(key, "a", Role::Contributor),
        credential(other, "b", Role::Reader),
    ])
    .unwrap();
    let mirror = registry.clone();
    assert_eq!(
        mirror
            .authenticate(header.as_bytes())
            .unwrap()
            .author
            .id
            .as_str(),
        "a"
    );
    assert_eq!(
        registry.authenticate(other_header.as_bytes()).unwrap().role,
        Role::Reader
    );
    assert!(registry.authenticate(b"Bearer invented").is_none());
    let a = AccessKey::read(&t.0.join("a.key")).unwrap();
    let same = AccessKey::read(&t.0.join("a.key")).unwrap();
    assert!(
        registry
            .replace(vec![
                credential(a, "a", Role::Contributor),
                credential(same, "other", Role::Maintainer)
            ])
            .is_err()
    );
    assert!(registry.authenticate(header.as_bytes()).is_some());
    let a = AccessKey::read(&t.0.join("a.key")).unwrap();
    let b = AccessKey::read(&t.0.join("b.key")).unwrap();
    assert!(
        registry
            .replace(vec![
                credential(a, "duplicate", Role::Contributor),
                credential(b, "duplicate", Role::Maintainer)
            ])
            .is_err()
    );
    let b = AccessKey::read(&t.0.join("b.key")).unwrap();
    registry
        .replace(vec![credential(b, "b", Role::Maintainer)])
        .unwrap();
    assert!(mirror.authenticate(header.as_bytes()).is_none());
    assert_eq!(
        mirror.authenticate(other_header.as_bytes()).unwrap().role,
        Role::Maintainer
    );
}
#[test]
fn reader_permissions_cover_all_current_document_mutations_and_reliable_apply() {
    let t = Temp::new();
    let (key, header) = t.key("reader.key");
    let registry = MemberRegistry::new(vec![credential(key, "r", Role::Reader)]).unwrap();
    let reader = registry.authenticate(header.as_bytes()).unwrap();
    let folder = wire::FolderId::new(format!("fld-{}", "a".repeat(32))).unwrap();
    let id = wire::DocumentId::new(format!("doc-{}", "b".repeat(32))).unwrap();
    let rev = wire::RevisionId::new("c".repeat(64)).unwrap();
    let name = wire::Text::new("name").unwrap();
    let payload = wire::Payload {
        kind: wire::PayloadKind::NativeDocument,
        bytes: wire::Blob::new(vec![1]).unwrap(),
        archive_hash: wire::Digest::of(&[1]),
        selection: None,
        sharing: None,
    };
    let mutations = vec![
        wire::Request::CreateFolder {
            parent: None,
            name: name.clone(),
        },
        wire::Request::RenameFolder {
            id: folder.clone(),
            expected_generation: 1,
            name: name.clone(),
        },
        wire::Request::MoveFolder {
            id: folder.clone(),
            expected_generation: 1,
            parent: None,
        },
        wire::Request::DeleteEmptyFolder {
            id: folder,
            expected_generation: 1,
        },
        wire::Request::ImportDocument {
            folder: None,
            name: name.clone(),
            payload: payload.clone(),
        },
        wire::Request::SaveRevision {
            id: id.clone(),
            base_revision: rev,
            payload,
        },
        wire::Request::RenameDocument {
            id: id.clone(),
            expected_generation: 1,
            name: name.clone(),
        },
        wire::Request::MoveDocument {
            id: id.clone(),
            expected_generation: 1,
            folder: None,
        },
    ];
    for request in mutations {
        assert!(!reader.permits_library(&request));
    }
    for request in [
        wire::Request::Browse {
            parent: None,
            page: Default::default(),
        },
        wire::Request::OpenDocument {
            id: id.clone(),
            revision: None,
        },
        wire::Request::Search {
            query: wire::Text::new("p").unwrap(),
            scope: wire::SearchScope::AllFolders,
            page: Default::default(),
        },
        wire::Request::ListRevisions {
            id,
            page: Default::default(),
        },
    ] {
        assert!(reader.permits_library(&request));
    }
    use collagene_library_sync as sync;
    let mutation = sync::Mutation {
        id: sync::MutationId::new(format!("mut-{}", "d".repeat(32))).unwrap(),
        operation: sync::MutationOperation::CreateFolder { parent: None, name },
    };
    assert!(!reader.permits_sync(&sync::Request::Apply {
        mutation: mutation.clone()
    }));
    assert!(reader.permits_sync(&sync::Request::Lookup {
        mutation_id: mutation.id.clone(),
        expected_intent: wire::Digest::of(b"intent")
    }));
    let (key, header) = t.key("contributor.key");
    registry
        .replace(vec![credential(key, "c", Role::Contributor)])
        .unwrap();
    assert!(
        registry
            .authenticate(header.as_bytes())
            .unwrap()
            .permits_sync(&sync::Request::Apply { mutation })
    );
}
#[test]
fn operator_registry_rejects_untrusted_fields_versions_exposed_secrets_and_unsafe_files() {
    let t = Temp::new();
    let assets = t.0.join("assets");
    std::fs::create_dir(&assets).unwrap();
    let (_, header) = t.key("a.key");
    let config = r#"{"schema_version":1,"members":[{"actor_id":"a","display_name":"Scientist","role":"Contributor","key_file":"a.key"}]}"#;
    let path = t.config(config);
    let registry = MemberRegistry::read_file(&path, &assets).unwrap();
    assert_eq!(
        registry
            .authenticate(header.as_bytes())
            .unwrap()
            .author
            .display_name,
        "Scientist"
    );
    for invalid in [
        config.replace("\"schema_version\":1", "\"schema_version\":2"),
        config.replace(
            "\"role\":\"Contributor\"",
            "\"role\":\"Contributor\",\"verified\":true",
        ),
        config.replace("\"actor_id\":\"a\"", "\"actor_id\":\"fake identity\""),
    ] {
        t.config(&invalid);
        assert!(MemberRegistry::read_file(&path, &assets).is_err());
    }
    t.config(config);
    assert!(MemberRegistry::read_file(&path, &t.0).is_err());
    std::fs::copy(t.0.join("a.key"), assets.join("public.key")).unwrap();
    t.config(&config.replace("a.key", "assets/public.key"));
    assert!(MemberRegistry::read_file(&path, &assets).is_err());
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        t.config(config);
        std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o644)).unwrap();
        assert!(MemberRegistry::read_file(&path, &assets).is_err());
    }
}
