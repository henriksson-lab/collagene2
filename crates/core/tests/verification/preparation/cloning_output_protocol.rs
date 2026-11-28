use collagene_core::{
    cloning::restriction::{
        catalog as c,
        cut_selection::{self as cuts, Boundary, Choice},
    },
    formats::oligos::Delimiter,
    planning::{
        oligos as plan_oligos,
        plan::{DesignPlan, Operation, PlanInput, Source},
    },
    protocol::{
        self,
        request::Request,
        schema::{AnalysisJob, AnalysisResult},
    },
    sequence::analysis::Primer,
};
use collagene_model::{Molecule, Topology};
fn wire(version: u16, request: &Request) -> Vec<u8> {
    let p = postcard::to_allocvec(request).unwrap();
    let mut b = b"CGQ\0".to_vec();
    b.extend(version.to_le_bytes());
    b.extend((p.len() as u32).to_le_bytes());
    b.extend(blake3::hash(&p).as_bytes());
    b.extend(p);
    b
}
fn response(request: &Request, tag: u32, result_tag: u32) -> AnalysisResult {
    let bytes = protocol::request::Request::encode(request).unwrap();
    let ((_, _, actual), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(actual, tag);
    assert!(protocol::execution::process_request(&wire(37, request)).is_err());
    let bytes = protocol::execution::process_request(&bytes).unwrap();
    let ((_, _, actual), _) = postcard::take_from_bytes::<(u32, &str, u32)>(&bytes[42..]).unwrap();
    assert_eq!(actual, result_tag);
    let mut old = bytes.clone();
    old[4..6].copy_from_slice(&37u16.to_le_bytes());
    assert!(protocol::response::Response::decode(&old).is_err());
    let decoded = protocol::response::Response::decode(&bytes).unwrap();
    assert_eq!(decoded.id, request.id);
    assert_eq!(decoded.revision, request.revision);
    decoded.result
}
#[test]
fn two_cut_choice_reconstructs_exact_fragment_and_binds_the_selected_path() {
    let catalog = c::builtin_catalog().unwrap();
    let selection = c::select(&catalog, &["ecori".into()]).unwrap();
    let molecule = Molecule::dna(
        "two actual sites",
        "GAATTCAAAAGAATTCCCC",
        Topology::Circular,
    )
    .unwrap();
    let options = collagene_core::cloning::restriction::inspection::InspectionOptions::default();
    let batch = c::batch_revision(
        &catalog,
        &selection,
        "source-event-role",
        &molecule,
        &options,
    )
    .unwrap();
    let choice = Choice {
        first: Boundary::Cut { site_index: 0 },
        second: Boundary::Cut { site_index: 1 },
        complement: false,
    };
    let revision = cuts::selection_revision(
        &catalog,
        &selection,
        "source-event-role",
        &molecule,
        &options,
        &batch,
        &choice,
    )
    .unwrap();
    let request = Request {
        id: 180,
        revision: revision.clone(),
        job: AnalysisJob::SelectDigestFragment {
            selection: selection.clone(),
            source_role: "source-event-role".into(),
            molecule: molecule.clone(),
            options: options.clone(),
            expected_batch_revision: batch.clone(),
            choice: choice.clone(),
        },
    };
    let AnalysisResult::DigestFragmentSelected(report) = response(&request, 97, 89) else {
        panic!("expected fresh cut selection")
    };
    assert_eq!(
        report.product.product.molecule.sequence.as_str(),
        "AATTCAAAAG"
    );
    report.validate_binding(&revision).unwrap();
    let changed = Request {
        id: 180,
        revision,
        job: AnalysisJob::SelectDigestFragment {
            selection,
            source_role: "source-event-role".into(),
            molecule,
            options,
            expected_batch_revision: batch,
            choice: Choice {
                complement: true,
                ..choice
            },
        },
    };
    assert!(
        protocol::execution::process_request(
            &protocol::request::Request::encode(&changed).unwrap()
        )
        .is_err()
    );
}
fn pcr() -> (DesignPlan, Molecule) {
    let molecule = Molecule::dna(
        "literal template",
        "AACCGGTTAGCTACTAGCTATGCCCCTTGGAA",
        Topology::Linear,
    )
    .unwrap();
    let target = Molecule::dna(
        "independent tailed PCR",
        "CGAACCGGTTAGCTACTAGCTATGCCCCTTGGAA",
        Topology::Linear,
    )
    .unwrap();
    let plan = DesignPlan {
        standard: None,
        name: "ordered primer graph".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "template".into(),
            revision: molecule.revision(),
            molecule,
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("template".into()),
            forward: Primer::new("=F", "AACCGGTT", "CG").unwrap(),
            reverse: Primer::new("R", "TTCCAAGG", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        }],
    };
    (plan, target)
}
#[test]
fn order_sheet_and_html_are_exported_from_the_bound_plan_and_target() {
    let (plan, target) = pcr();
    let options = plan_oligos::Options::default();
    let revision =
        plan_oligos::export_revision(&plan, &target, &options, Delimiter::Comma).unwrap();
    let request = Request {
        id: 181,
        revision,
        job: AnalysisJob::ExportPlanOligos {
            plan: plan.clone(),
            target: target.clone(),
            options,
            delimiter: Delimiter::Comma,
        },
    };
    let AnalysisResult::PlanOligosExported(sheet) = response(&request, 98, 90) else {
        panic!("expected oligo sheet")
    };
    assert!(sheet.report.complete);
    assert_eq!(
        sheet.report.plan_revision,
        collagene_core::planning::plan::DesignPlan::revision(&plan)
    );
    assert_eq!(sheet.report.oligos.records.len(), 2);
    assert!(sheet.export.text.contains("CGAACCGGTT"));
    assert!(sheet.export.text.contains("'=F"));
    assert_eq!(sheet.report.uses[0].primer.as_ref().unwrap().tail, "CG");
    let request = Request {
        id: 182,
        revision: target.revision(),
        job: AnalysisJob::ExportProtocolHtml {
            plan: plan.clone(),
            target: target.clone(),
        },
    };
    let AnalysisResult::ProtocolHtml(html) = response(&request, 99, 91) else {
        panic!("expected HTML steps")
    };
    assert!(html.contains("CGAACCGGTT"));
    assert!(html.contains(&collagene_core::planning::plan::DesignPlan::revision(&plan)));
}
