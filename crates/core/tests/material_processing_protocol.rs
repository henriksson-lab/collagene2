//! Current binary workers preserve target-free mismatch PCR and actual-strand repair.
use collagene_core::{
    Error, Molecule, Topology,
    analysis::Primer,
    catalog_digest,
    chemistry::EndRepairMethod,
    material_derivation as native,
    mutagenesis::PcrOptions,
    protocol::{self, AnalysisJob, AnalysisResult, Request},
    run_lineage::EventBinding,
    run_material::MaterialOrigin,
    workspace::{self, Workspace},
};
const STOCK: &str = "AGTCCGATGAATTCCGTAGGACTTCAGACCTGTAGATCGTCA";
const MUTANT: &str = "CCAGTTCGATGAATTCCGTAGGACTTCAGACCTGTAGATCGTCACC";
fn event(w: &Workspace) -> EventBinding {
    let d = &w.documents[0];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn wire(magic: &[u8; 4], payload: &impl serde::Serialize) -> Vec<u8> {
    let p = postcard::to_allocvec(payload).unwrap();
    let mut b = magic.to_vec();
    b.extend(protocol::WORKER_SCHEMA_VERSION.to_le_bytes());
    b.extend((p.len() as u32).to_le_bytes());
    b.extend(blake3::hash(&p).as_bytes());
    b.extend(p);
    b
}
fn tag(v: &impl serde::Serialize) -> u32 {
    postcard::take_from_bytes::<u32>(&postcard::to_allocvec(v).unwrap())
        .unwrap()
        .0
}
fn call(job: AnalysisJob, revision: String, selector: u32, result: u32) -> protocol::Response {
    assert_eq!(tag(&job), selector);
    let q = Request {
        id: selector,
        revision: revision.clone(),
        job,
    };
    let b = protocol::encode_request(&q).unwrap();
    assert_eq!(&b[4..6], &protocol::WORKER_SCHEMA_VERSION.to_le_bytes());
    let b = protocol::AnalysisSession::default()
        .process_request(&b)
        .unwrap();
    assert_eq!(&b[4..6], &protocol::WORKER_SCHEMA_VERSION.to_le_bytes());
    let r = protocol::decode_response(&b).unwrap();
    assert_eq!(r.revision, revision);
    assert_eq!(tag(&r.result), result, "{:?}", r.result);
    r
}
#[test]
fn mismatch_install_redigest_and_repair_keep_literal_tails_and_actual_sticky_strands() {
    let source = Molecule::dna("Literal original", STOCK, Topology::Linear).unwrap();
    let w = Workspace::from_documents("Original", vec![source.clone()], 0).unwrap();
    let before = workspace::save_workspace(&w).unwrap();
    let forward = Primer::new("F mismatch", "AGTTCGAT", "CC").unwrap();
    let reverse = Primer::new("R exact", "TGACGATC", "GG").unwrap();
    let mut options = PcrOptions::default();
    options.binding.min_three_prime_exact = 4;
    let origin = native::SourceOrigin::OriginalReference;
    let context =
        native::mismatch_pcr_preview_revision(&source, &origin, &forward, &reverse, &options)
            .unwrap();
    let job = AnalysisJob::InspectMaterialMismatchPcr {
        source: source.clone(),
        origin: origin.clone(),
        forward: forward.clone(),
        reverse: reverse.clone(),
        options: options.clone(),
    };
    let response = call(job.clone(), context.clone(), 175, 165);
    let AnalysisResult::MaterialMismatchPcrPreview(preview) = response.result else {
        panic!("mismatch preview")
    };
    preview.validate_binding(&context).unwrap();
    assert!(preview.reusable());
    assert_eq!(preview.preview.documents[0].sequence.as_str(), MUTANT);
    assert_eq!(preview.preview.report.products[0].substitutions.len(), 1);
    assert_eq!(
        preview.preview.report.products[0].substitutions[0].source_position,
        3
    );
    let mut changed = job;
    let AnalysisJob::InspectMaterialMismatchPcr {
        forward: changed_forward,
        ..
    } = &mut changed
    else {
        unreachable!()
    };
    *changed_forward = Primer::new("F mismatch", "AGTTCGAT", "AA").unwrap();
    assert!(
        protocol::process_request(
            &protocol::encode_request(&Request {
                id: 1,
                revision: context,
                job: changed
            })
            .unwrap()
        )
        .is_err()
    );
    let request = native::Request {
        sources: vec![collagene_core::material_derivation::BoundSource {
            role: collagene_core::material_derivation::SourceRole::Template,
            event: event(&w),
            origin: origin,
        }],
        derivation: native::Derivation::MismatchPcr {
            forward,
            reverse,
            options,
            selection: native::MismatchPcrSelection {
                selected_index: 0,
                expected_preview_revision: preview.preview.documents[0].revision(),
                geometry: preview.preview.report.products[0].clone(),
            },
        },
        name: "Reusable deliberate substitution".into(),
        save_primers: true,
    };
    let expected = native::revision(&w, &request).unwrap();
    let response = call(
        AnalysisJob::InstallMaterialDerivation {
            workspace: w.clone(),
            request,
        },
        expected,
        171,
        161,
    );
    let AnalysisResult::MaterialDerivationInstalled(installed) = response.result else {
        panic!("installation")
    };
    let pcr = workspace::load_workspace(&installed.archive).unwrap();
    installed.receipt.validate_archive(&pcr).unwrap();
    assert!(installed.receipt.primer_attachment.is_some());
    assert_eq!(
        native::reconstruct(&pcr, &installed.receipt.product)
            .unwrap()
            .reference()
            .sequence
            .as_str(),
        MUTANT
    );
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let eco = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap();
    let request = native::Request {
        sources: vec![collagene_core::material_derivation::BoundSource {
            role: collagene_core::material_derivation::SourceRole::Template,
            event: installed.receipt.product.clone(),
            origin: native::SourceOrigin::Saved(MaterialOrigin::Derived(
                native::capsule(&pcr, &installed.receipt.product).unwrap(),
            )),
        }],
        derivation: native::Derivation::CatalogDigest {
            selection: catalog_digest::select(&catalog, &[eco.id.clone()]).unwrap(),
            options: Default::default(),
            fragment_index: 1,
        },
        name: "Actual EcoRI right fragment".into(),
        save_primers: false,
    };
    let cut = native::prepare_install(&pcr, &request).unwrap();
    let source = cut.workspace.revisions[&cut.receipt.product.revision].clone();
    let origin = native::SourceOrigin::Saved(MaterialOrigin::Derived(
        native::capsule(&cut.workspace, &cut.receipt.product).unwrap(),
    ));
    let context =
        native::repair_preview_revision(&source, &origin, EndRepairMethod::FillFivePrimeOverhangs)
            .unwrap();
    let response = call(
        AnalysisJob::InspectMaterialRepair {
            source: source.clone(),
            origin: origin.clone(),
            method: EndRepairMethod::FillFivePrimeOverhangs,
        },
        context.clone(),
        176,
        166,
    );
    let AnalysisResult::MaterialRepairPreview(repair) = response.result else {
        panic!("repair")
    };
    repair.validate_binding(&context).unwrap();
    assert_eq!(repair.repair.source.bottom_start, 4);
    assert_eq!(repair.repair.duplex.bottom_start, 0);
    assert_eq!(repair.repair.bottom_added_3prime, "AATT");
    assert!(repair.repair.mappings.iter().any(|m| m.synthesized));
    let mut forged = protocol::Response {
        id: 181,
        revision: context.clone(),
        result: AnalysisResult::MaterialRepairPreview(repair),
    };
    let AnalysisResult::MaterialRepairPreview(report) = &mut forged.result else {
        unreachable!()
    };
    report.context_revision = "0".repeat(64);
    assert!(protocol::decode_response(&wire(b"CGR\0", &forged)).is_err());
    let request = native::Request {
        sources: vec![collagene_core::material_derivation::BoundSource {
            role: collagene_core::material_derivation::SourceRole::Template,
            event: cut.receipt.product.clone(),
            origin: origin,
        }],
        derivation: native::Derivation::EndRepair {
            method: EndRepairMethod::FillFivePrimeOverhangs,
        },
        name: "Filled saved fragment".into(),
        save_primers: false,
    };
    let expected = native::revision(&cut.workspace, &request).unwrap();
    let response = call(
        AnalysisJob::InstallMaterialDerivation {
            workspace: cut.workspace,
            request,
        },
        expected,
        171,
        161,
    );
    let AnalysisResult::MaterialDerivationInstalled(installed) = response.result else {
        panic!("repaired installation")
    };
    let repaired = workspace::load_workspace(&installed.archive).unwrap();
    let fresh = native::reconstruct(&repaired, &installed.receipt.product).unwrap();
    assert_eq!(fresh.tracked().duplex.bottom_start, 0);
    assert_eq!(
        fresh.reference().sequence.as_str(),
        "AATTCCGTAGGACTTCAGACCTGTAGATCGTCACC"
    );
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
}
#[test]
fn preview_admission_uses_compound_budget_but_map_remains_ordinary() {
    let mut body = postcard::to_allocvec(&(1u32, "context", 175u32)).unwrap();
    body.resize(collagene_core::formats::MAX_IMPORT_BYTES + 1, 0);
    for selector in [175u32, 176, 174] {
        let prefix = postcard::to_allocvec(&(1u32, "context", selector)).unwrap();
        body[..prefix.len()].copy_from_slice(&prefix);
        let mut b = b"CGQ\0".to_vec();
        b.extend(protocol::WORKER_SCHEMA_VERSION.to_le_bytes());
        b.extend((body.len() as u32).to_le_bytes());
        b.extend(blake3::hash(&body).as_bytes());
        b.extend(&body);
        let stateless = protocol::process_request(&b).unwrap_err();
        let session = protocol::AnalysisSession::default()
            .process_request(&b)
            .unwrap_err();
        for error in [stateless, session] {
            if selector == 174 {
                assert!(matches!(error, Error::TooLarge), "{error}");
            } else {
                assert!(
                    matches!(error, Error::Format(_)),
                    "compound preview was refused before typed decoding: {error}"
                );
            }
        }
    }
}
