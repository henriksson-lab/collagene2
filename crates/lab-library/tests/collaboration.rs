use collagene_lab_library::{collaboration::*, *};
use std::{
    path::PathBuf,
    sync::{Arc, Barrier},
};
fn payload(text: &str) -> Payload {
    Payload {
        kind: PayloadKind::NativeDocument,
        bytes: text.as_bytes().to_vec(),
        archive_hash: Digest::parse(blake3::hash(text.as_bytes()).to_hex().to_string()).unwrap(),
        selection: None,
        sharing: None,
    }
}
fn member(db: &Library, lab: &LabAccess, id: &str, role: Role) -> MemberAccess {
    db.trusted_member_access(
        lab,
        Author::new(ActorId::parse(id).unwrap(), format!("Scientist {id}")).unwrap(),
        role,
    )
    .unwrap()
}
fn fixture() -> (Library, LabAccess, MemberAccess, MemberAccess, Target) {
    let mut db = Library::open_in_memory().unwrap();
    let lab = db.create_lab("A").unwrap();
    let lab = db.trusted_host_access(&lab.id).unwrap();
    let a = member(&db, &lab, "a", Role::Contributor);
    let b = member(&db, &lab, "b", Role::Contributor);
    let d = db
        .create_document(&lab, None, "pExample", &payload("first"))
        .unwrap();
    let target = Target {
        document: d.id,
        revision: d.head.id,
    };
    (db, lab, a, b, target)
}
fn open(body: &str) -> Action {
    Action::OpenThread { body: body.into() }
}
#[test]
fn human_review_is_exact_revision_bound_and_never_changes_documents_or_catalog() {
    let (mut db, lab, a, b, target) = fixture();
    let before = db.lab(&lab).unwrap();
    let original = db.document(&lab, &target.document).unwrap();
    let first = db
        .append_discussion(
            &a,
            &target,
            0,
            Action::Review {
                decision: ReviewDecision::Approved,
                body: "Checked cloning design".into(),
            },
        )
        .unwrap();
    assert_eq!(first.author, *a.author());
    assert_eq!(db.document(&lab, &target.document).unwrap(), original);
    assert_eq!(db.lab(&lab).unwrap(), before);
    db.append_discussion(
        &b,
        &target,
        1,
        Action::Review {
            decision: ReviewDecision::ChangesRequested,
            body: "Check orientation".into(),
        },
    )
    .unwrap();
    db.append_discussion(
        &a,
        &target,
        2,
        Action::Review {
            decision: ReviewDecision::Withdrawn,
            body: String::new(),
        },
    )
    .unwrap();
    let new = db
        .save(&lab, &target.document, &target.revision, &payload("second"))
        .unwrap();
    let new_target = Target {
        document: new.id,
        revision: new.head.id,
    };
    let page = db.discussion_page(&a, &new_target, 0, 200, None).unwrap();
    assert_eq!(page.generation, 0);
    assert!(page.events.is_empty());
    let old = db.discussion_page(&b, &target, 0, 200, None).unwrap();
    assert_eq!(old.events.len(), 3);
    assert_eq!(old.events[0], first);
    assert_eq!(
        db.payload(&lab, &target.document, &target.revision)
            .unwrap()
            .bytes,
        b"first"
    );
}
#[test]
fn threads_enforce_authorship_resolution_and_immutable_history() {
    let (mut db, lab, a, b, target) = fixture();
    let admin = member(&db, &lab, "moderator", Role::Maintainer);
    let reader = member(&db, &lab, "reader", Role::Reader);
    db.append_discussion(&a, &target, 0, open("Verify the junction"))
        .unwrap();
    db.append_discussion(
        &b,
        &target,
        1,
        Action::Reply {
            thread: 1,
            body: "Literal bases agree".into(),
        },
    )
    .unwrap();
    assert!(matches!(
        db.append_discussion(&reader, &target, 2, open("x")),
        Err(DiscussionError::Forbidden)
    ));
    assert!(matches!(
        db.append_discussion(
            &b,
            &target,
            2,
            Action::SetResolved {
                thread: 1,
                resolved: true
            }
        ),
        Err(DiscussionError::Forbidden)
    ));
    assert!(
        db.append_discussion(
            &a,
            &target,
            2,
            Action::Reply {
                thread: 2,
                body: "Reply is not a root".into()
            }
        )
        .is_err()
    );
    db.append_discussion(
        &a,
        &target,
        2,
        Action::SetResolved {
            thread: 1,
            resolved: true,
        },
    )
    .unwrap();
    assert!(
        db.append_discussion(
            &b,
            &target,
            3,
            Action::Reply {
                thread: 1,
                body: "Cannot reply while resolved".into()
            }
        )
        .is_err()
    );
    assert!(
        db.append_discussion(
            &a,
            &target,
            3,
            Action::SetResolved {
                thread: 1,
                resolved: true
            }
        )
        .is_err()
    );
    db.append_discussion(
        &admin,
        &target,
        3,
        Action::SetResolved {
            thread: 1,
            resolved: false,
        },
    )
    .unwrap();
    db.append_discussion(
        &b,
        &target,
        4,
        Action::Reply {
            thread: 1,
            body: "Now reopened".into(),
        },
    )
    .unwrap();
    let events = db
        .discussion_page(&reader, &target, 0, 200, None)
        .unwrap()
        .events;
    assert_eq!(events.len(), 5);
    assert_eq!(events[0].author.id.as_str(), "a");
    assert_eq!(events[3].role, Role::Maintainer);
}
#[test]
fn scoped_revision_lookup_snapshot_pagination_and_lost_ack_refuse_atomically() {
    let (mut db, lab, a, _, target) = fixture();
    let other = db.create_lab("B").unwrap();
    let other = db.trusted_host_access(&other.id).unwrap();
    let foreign = member(&db, &other, "a", Role::Maintainer);
    assert!(matches!(
        db.discussion_page(&foreign, &target, 0, 10, None),
        Err(DiscussionError::Storage(Error::NotFound))
    ));
    assert!(
        db.append_discussion(&foreign, &target, 0, open("foreign"))
            .is_err()
    );
    let wrong = db
        .create_document(&lab, None, "Other", &payload("other"))
        .unwrap();
    let wrong = Target {
        document: target.document.clone(),
        revision: wrong.head.id,
    };
    assert!(
        db.append_discussion(&a, &wrong, 0, open("wrong revision"))
            .is_err()
    );
    db.append_discussion(&a, &target, 0, open("first")).unwrap();
    assert!(matches!(
        db.append_discussion(&a, &target, 0, open("first")),
        Err(DiscussionError::Conflict {
            expected: 0,
            actual: 1
        })
    ));
    db.append_discussion(&a, &target, 1, open("second"))
        .unwrap();
    let first = db.discussion_page(&a, &target, 0, 1, None).unwrap();
    assert_eq!(first.generation, 2);
    assert_eq!(first.next_offset, Some(1));
    assert!(db.discussion_page(&a, &target, 1, 1, None).is_err());
    assert_eq!(
        db.discussion_page(&a, &target, 1, 1, Some(2))
            .unwrap()
            .events[0]
            .sequence,
        2
    );
    db.append_discussion(&a, &target, 2, open("third")).unwrap();
    assert!(matches!(
        db.discussion_page(&a, &target, 1, 1, Some(2)),
        Err(DiscussionError::Conflict {
            expected: 2,
            actual: 3
        })
    ));
    assert_eq!(
        db.discussion_page(&a, &target, 0, 10, None)
            .unwrap()
            .events
            .len(),
        3
    );
}
#[test]
fn input_bounds_are_bytes_and_author_claims_are_not_taken_from_comment_text() {
    let (mut db, lab, a, _, target) = fixture();
    for body in [
        " ".into(),
        "bad\0text".into(),
        "é".repeat(MAX_BODY_BYTES / 2 + 1),
    ] {
        assert!(
            db.append_discussion(&a, &target, 0, Action::OpenThread { body })
                .is_err()
        );
    }
    assert!(ActorId::parse("not an identity").is_err());
    assert!(Author::new(ActorId::parse("a").unwrap(), "bad\nname").is_err());
    let mut invalid = a.author().clone();
    invalid.display_name = String::new();
    assert!(
        db.trusted_member_access(&lab, invalid, Role::Maintainer)
            .is_err()
    );
    let body = "<script>author=moderator; approved=true</script>\n\t";
    let e = db.append_discussion(&a, &target, 0, open(body)).unwrap();
    assert_eq!(e.author.id.as_str(), "a");
    assert_eq!(e.role, Role::Contributor);
    assert_eq!(e.action, open(body)); // Plain text; UI must escape, not interpret.
    db.append_discussion(&a, &target, 1, open(&"é".repeat(MAX_BODY_BYTES / 2)))
        .unwrap();
}
struct Temp(PathBuf);
impl Temp {
    fn new() -> Self {
        let p = std::env::temp_dir().join(format!(
            "collagene-discussion-{}-{}",
            std::process::id(),
            std::time::SystemTime::now()
                .duration_since(std::time::UNIX_EPOCH)
                .unwrap()
                .as_nanos()
        ));
        let mut builder = std::fs::DirBuilder::new();
        #[cfg(unix)]
        {
            use std::os::unix::fs::DirBuilderExt;
            builder.mode(0o700);
        }
        builder.create(&p).unwrap();
        Self(p)
    }
    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}
impl Drop for Temp {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}
#[test]
fn independent_connections_admit_one_concurrent_writer_and_retain_events_on_restart() {
    let temp = Temp::new();
    let path = temp.path("db.sqlite");
    let mut db = Library::open(&path).unwrap();
    let lab = db.create_lab("A").unwrap();
    let access = db.trusted_host_access(&lab.id).unwrap();
    let doc = db
        .create_document(&access, None, "pExample", &payload("first"))
        .unwrap();
    let target = Target {
        document: doc.id,
        revision: doc.head.id,
    };
    let a = member(&db, &access, "a", Role::Contributor);
    let barrier = Arc::new(Barrier::new(2));
    let handles: Vec<_> = ["a", "b"]
        .into_iter()
        .map(|id| {
            let barrier = barrier.clone();
            let path = path.clone();
            let access = access.clone();
            let target = target.clone();
            std::thread::spawn(move || {
                let mut db = Library::open(path).unwrap();
                let m = member(&db, &access, id, Role::Contributor);
                barrier.wait();
                db.append_discussion(&m, &target, 0, open(id))
            })
        })
        .collect();
    let outcomes: Vec<_> = handles.into_iter().map(|h| h.join().unwrap()).collect();
    assert_eq!(outcomes.iter().filter(|r| r.is_ok()).count(), 1);
    assert_eq!(
        outcomes
            .iter()
            .filter(|r| matches!(
                r,
                Err(DiscussionError::Conflict {
                    expected: 0,
                    actual: 1
                })
            ))
            .count(),
        1
    );
    drop(db);
    let mut db = Library::open(&path).unwrap();
    assert_eq!(
        db.discussion_page(&a, &target, 0, 10, None)
            .unwrap()
            .events
            .len(),
        1
    );
    let raw = rusqlite::Connection::open(&path).unwrap();
    assert!(
        raw.execute("UPDATE discussion_events SET actor_id='forged'", [])
            .is_err()
    );
    assert!(raw.execute("DELETE FROM discussion_events", []).is_err());
    // Quota accounting is bounded independently from sequence archive size.
    raw.execute("UPDATE labs SET discussion_bytes=16777216", [])
        .unwrap();
    assert!(matches!(
        db.append_discussion(&a, &target, 1, open("over quota")),
        Err(DiscussionError::Storage(Error::Limit(
            "lab discussion quota"
        )))
    ));
    assert_eq!(
        db.discussion_page(&a, &target, 0, 10, None)
            .unwrap()
            .generation,
        1
    );
    let policy = backup::Policy::single_lab(lab.id);
    assert!(backup::create(&db, &temp.path("corrupt.cglbak"), &policy).is_err());
}
#[test]
fn backup_restores_discussion_and_rejects_superseded_schema_without_modifying_it() {
    let temp = Temp::new();
    let path = temp.path("db.sqlite");
    let mut db = Library::open(&path).unwrap();
    let lab = db.create_lab("A").unwrap();
    let access = db.trusted_host_access(&lab.id).unwrap();
    let a = member(&db, &access, "a", Role::Contributor);
    let doc = db
        .create_document(&access, None, "pExample", &payload("first"))
        .unwrap();
    let target = Target {
        document: doc.id,
        revision: doc.head.id,
    };
    db.append_discussion(&a, &target, 0, open("Junction review"))
        .unwrap();
    db.append_discussion(
        &a,
        &target,
        1,
        Action::SetResolved {
            thread: 1,
            resolved: true,
        },
    )
    .unwrap();
    db.append_discussion(
        &a,
        &target,
        2,
        Action::SetResolved {
            thread: 1,
            resolved: false,
        },
    )
    .unwrap();
    db.append_discussion(
        &a,
        &target,
        3,
        Action::Reply {
            thread: 1,
            body: "Additional check".into(),
        },
    )
    .unwrap();
    db.append_discussion(
        &a,
        &target,
        4,
        Action::Review {
            decision: ReviewDecision::Approved,
            body: String::new(),
        },
    )
    .unwrap();
    let original = db.discussion_page(&a, &target, 0, 200, None).unwrap();
    let policy = backup::Policy::single_lab(lab.id);
    let backup_file = temp.path("backup.cglbak");
    let manifest = backup::create(&db, &backup_file, &policy).unwrap();
    assert_eq!(manifest.schema_version, 2);
    assert_eq!(manifest.table_rows["discussion_events"], 5);
    backup::restore_new(&backup_file, &temp.path("restored"), &policy).unwrap();
    let mut restored = Library::open(temp.path("restored/library.sqlite")).unwrap();
    assert_eq!(
        restored.discussion_page(&a, &target, 0, 200, None).unwrap(),
        original
    );
    drop(db);
    let raw = rusqlite::Connection::open(&path).unwrap();
    raw.execute_batch("PRAGMA user_version=1;").unwrap();
    drop(raw);
    assert!(matches!(
        Library::open(&path),
        Err(Error::UnsupportedDatabase)
    ));
    let raw = rusqlite::Connection::open(&path).unwrap();
    assert_eq!(
        raw.query_row("PRAGMA user_version", [], |r| r.get::<_, i64>(0))
            .unwrap(),
        1
    );
    assert_eq!(
        raw.query_row("SELECT count(*) FROM discussion_events", [], |r| r
            .get::<_, u64>(0))
            .unwrap(),
        5
    );
}
