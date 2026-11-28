use collagene_core::{
    history::DocumentHistory,
    oligo_formats::{self, QualityEncoding},
    read_checks,
    read_review::*,
    run_lineage::EventBinding,
    trace_validation::{TraceValidationOptions, TrimMode},
    workspace::Workspace,
};
use collagene_model::{Molecule, Topology};
fn binding(w: &Workspace, n: usize) -> EventBinding {
    let d = &w.documents[n];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn fixture(reads: &[&str]) -> (Workspace, String) {
    fixture_at("TTTCGATACGGG", Topology::Linear, reads)
}
fn fixture_at(reference: &str, topology: Topology, reads: &[&str]) -> (Workspace, String) {
    let mut molecules = vec![Molecule::dna("construct", reference, topology).unwrap()];
    for (i, r) in reads.iter().enumerate() {
        let fastq = format!("@read{i}\n{r}\n+\n{}\n", "I".repeat(r.len()));
        molecules.push(
            oligo_formats::import_fastq(fastq.as_bytes(), QualityEncoding::Phred33)
                .unwrap()
                .molecules
                .remove(0),
        );
    }
    let w = Workspace::from_documents("explicit review", molecules, 0).unwrap();
    let mut options = TraceValidationOptions::default();
    options.trim.mode = TrimMode::AllCalls;
    options.trim.min_retained_bases = 1;
    let c = read_checks::capture(
        &w,
        binding(&w, 0),
        (1..w.documents.len()).map(|i| binding(&w, i)).collect(),
        &options,
    )
    .unwrap();
    (c.workspace, c.attachment_hash)
}
fn request(hash: &str, choice: Choice) -> Request {
    Request {
        attachment_hash: hash.into(),
        name: "Caller reviewed reference".into(),
        decisions: vec![Decision {
            key: EvidenceKey::Site { position: 5 },
            choice,
            reason: "Explicit review of the original eligible call".into(),
        }],
    }
}
fn allele() -> Choice {
    Choice::UseBase {
        read_index: 0,
        call_index: 5,
        base: 'T',
    }
}
#[test]
fn eligible_literal_allele_creates_new_guarded_reference_with_all_original_parents() {
    let (w, h) = fixture(&["TTTCGTTACGGG"]);
    let before = serde_json::to_value(&w).unwrap();
    let r = request(&h, allele());
    let p = preview(&w, &r).unwrap();
    p.validate_binding(&preview_revision(&w, &r).unwrap())
        .unwrap();
    assert_eq!(p.target.as_ref().unwrap().sequence.as_str(), "TTTCGTTACGGG");
    assert_eq!(p.edits[0].start, 5);
    assert_eq!(p.edits[0].replacement, "T");
    assert!(collagene_core::chemistry::Duplex::from_molecule(p.target.as_ref().unwrap()).is_err());
    let c = commit(&w, &r).unwrap();
    c.validate_binding(&commit_revision(&w, &r).unwrap())
        .unwrap();
    assert_eq!(serde_json::to_value(&w).unwrap(), before);
    let target = c.target.as_ref().unwrap();
    assert_eq!(c.workspace.documents.len(), w.documents.len() + 1);
    assert_eq!(
        c.workspace
            .head(&w.documents[0].id)
            .unwrap()
            .sequence
            .as_str(),
        "TTTCGATACGGG"
    );
    let event = c
        .workspace
        .history
        .iter()
        .find(|e| e.id == target.event_id)
        .unwrap();
    assert_eq!(
        event
            .parents
            .iter()
            .collect::<std::collections::BTreeSet<_>>(),
        [
            &w.documents[0].head_event.clone().unwrap(),
            &w.documents[1].head_event.clone().unwrap()
        ]
        .into_iter()
        .collect()
    );
    let saved = load(&c.workspace, &c.attachment_hash).unwrap();
    assert_eq!(saved.request, r);
    assert_eq!(saved.report.reference_revision, w.documents[0].head);
    assert_eq!(saved.target, c.target);
    assert_eq!(list(&c.workspace).unwrap().len(), 1);
    let again = commit(&c.workspace, &r).unwrap();
    assert_eq!(again.target, c.target);
    assert_eq!(again.attachment_hash, c.attachment_hash);
    assert_eq!(again.workspace.history.len(), c.workspace.history.len());
}
#[test]
fn retain_and_unresolved_record_decisions_without_changing_a_reference() {
    let (w, h) = fixture(&["TTTCGTTACGGG"]);
    let mut r = request(&h, Choice::RetainReference);
    r.decisions.push(Decision {
        key: EvidenceKey::Site { position: 0 },
        choice: Choice::Unresolved,
        reason: "No correction authorized".into(),
    });
    let p = preview(&w, &r).unwrap();
    assert!(p.target.is_none());
    assert!(p.edits.is_empty());
    let c = commit(&w, &r).unwrap();
    assert!(c.target.is_none());
    assert_eq!(c.workspace.documents.len(), w.documents.len());
    assert_eq!(c.workspace.history.len(), w.history.len());
    let again = commit(&c.workspace, &r).unwrap();
    assert_eq!(again.attachment_hash, c.attachment_hash);
}
#[test]
fn explicit_conflicting_allele_keeps_all_observations_and_reports_that_it_is_not_consensus() {
    let (w, h) = fixture(&["TTTCGTTACGGG", "TTTCGATACGGG"]);
    let r = request(&h, allele());
    let p = preview(&w, &r).unwrap();
    assert_eq!(p.target.as_ref().unwrap().sequence.as_str(), "TTTCGTTACGGG");
    assert!(p.warnings.iter().any(|w| w.contains("conflicting allele")));
    let c = commit(&w, &r).unwrap();
    let s = load(&c.workspace, &c.attachment_hash).unwrap();
    assert_eq!(s.report.reads.len(), 2);
    let t = c.target.unwrap();
    assert_eq!(
        c.workspace
            .history
            .iter()
            .find(|e| e.id == t.event_id)
            .unwrap()
            .parents
            .len(),
        3
    );
    let wrong = request(
        &h,
        Choice::UseBase {
            read_index: 0,
            call_index: 5,
            base: 'A',
        },
    );
    assert!(preview(&w, &wrong).is_err());
    let wrong = request(
        &h,
        Choice::UseBase {
            read_index: 0,
            call_index: 6,
            base: 'T',
        },
    );
    assert!(preview(&w, &wrong).is_err());
}
#[test]
fn eligible_insertions_and_deletions_use_original_coordinates_and_expected_full_calls() {
    let (w, h) = fixture(&["TTTCGACTACGGG", "TTTCGAGACGGG"]);
    let check = read_checks::load(&w, &h).unwrap();
    let insertion = &check.report.insertions[0];
    assert_eq!(insertion.sequence, "C");
    assert_eq!(insertion.reference_boundary, 6);
    let r = Request {
        attachment_hash: h.clone(),
        name: "Insertion and base choice".into(),
        decisions: vec![
            Decision {
                key: EvidenceKey::Insertion { index: 0 },
                choice: Choice::UseInsertion {
                    read_index: 0,
                    boundary: 6,
                    sequence: "C".into(),
                },
                reason: "Accept the complete observed insertion".into(),
            },
            Decision {
                key: EvidenceKey::Site { position: 6 },
                choice: Choice::UseBase {
                    read_index: 1,
                    call_index: 6,
                    base: 'G',
                },
                reason: "Choose this eligible conflicting call".into(),
            },
        ],
    };
    let p = preview(&w, &r).unwrap();
    assert_eq!(
        p.target.as_ref().unwrap().sequence.as_str(),
        "TTTCGACGACGGG"
    );
    assert_eq!((p.edits[0].start, p.edits[0].end), (6, 7));
    assert_eq!((p.edits[1].start, p.edits[1].end), (6, 6));
    let mut changed = r.clone();
    if let Choice::UseInsertion { sequence, .. } = &mut changed.decisions[0].choice {
        *sequence = "G".into();
    }
    assert!(preview(&w, &changed).is_err());
    let mut changed = r;
    if let Choice::UseInsertion { boundary, .. } = &mut changed.decisions[0].choice {
        *boundary = 5;
    }
    assert!(preview(&w, &changed).is_err());
    let (w, h) = fixture(&["TTTCGTACGGG"]);
    let p = preview(&w, &request(&h, Choice::DeleteBase { read_index: 0 })).unwrap();
    assert_eq!(p.target.as_ref().unwrap().sequence.as_str(), "TTTCGTACGGG");
    assert_eq!(p.edits[0].replacement, "");
    let (w, h) = fixture_at(
        "AGTCCGATGCTACGTAGGACTTCA",
        Topology::Circular,
        &["TTCACAGTC", "TTCAGAGTC"],
    );
    let check = read_checks::load(&w, &h).unwrap();
    assert_eq!(check.report.insertions.len(), 2);
    assert_eq!(check.report.insertions[0].reference_boundary % 24, 0);
    assert_eq!(check.report.insertions[1].reference_boundary % 24, 0);
    let decisions = check
        .report
        .insertions
        .iter()
        .enumerate()
        .map(|(index, i)| Decision {
            key: EvidenceKey::Insertion { index },
            choice: Choice::UseInsertion {
                read_index: i.read_index,
                boundary: i.reference_boundary,
                sequence: i.sequence.clone(),
            },
            reason: "Explicit competing insertion choice".into(),
        })
        .collect();
    let conflict = Request {
        attachment_hash: h,
        name: "Contradictory seam decisions".into(),
        decisions,
    };
    assert!(preview(&w, &conflict).is_err());
}
#[test]
fn edited_reviewed_output_is_preserved_and_new_derivation_is_idempotently_reused() {
    let (w, h) = fixture(&["TTTCGTTACGGG"]);
    let r = request(&h, allele());
    let first = commit(&w, &r).unwrap();
    let original = first.target.clone().unwrap();
    let mut changed = first.workspace.clone();
    let mut histories = DocumentHistory::restore_all(&changed).unwrap();
    histories
        .iter_mut()
        .find(|d| d.id == original.document_id)
        .unwrap()
        .commit(
            Molecule::dna("User edited reviewed output", "AAAA", Topology::Linear).unwrap(),
            "User edit after review",
        )
        .unwrap();
    changed.documents.clear();
    changed.history.clear();
    for d in &histories {
        d.contribute(&mut changed).unwrap();
    }
    let next = commit(&changed, &r).unwrap();
    assert_ne!(
        next.target.as_ref().unwrap().document_id,
        original.document_id
    );
    assert_eq!(
        next.workspace
            .head(&original.document_id)
            .unwrap()
            .sequence
            .as_str(),
        "AAAA"
    );
    assert_eq!(
        load(&next.workspace, &first.attachment_hash)
            .unwrap()
            .target,
        Some(original)
    );
    let again = commit(&next.workspace, &r).unwrap();
    assert_eq!(again.target, next.target);
    assert_eq!(again.attachment_hash, next.attachment_hash);
    assert_eq!(again.workspace.history.len(), next.workspace.history.len());
    let mut duplicate = r.clone();
    duplicate.decisions.push(duplicate.decisions[0].clone());
    assert!(source_preflight(&w, &duplicate).is_err());
    let mut reason = r.clone();
    reason.decisions[0].reason = " ".into();
    assert!(source_preflight(&w, &reason).is_err());
    let mut bounds = r;
    bounds.decisions[0].key = EvidenceKey::Site {
        position: usize::MAX,
    };
    assert!(source_preflight(&w, &bounds).is_err());
}
