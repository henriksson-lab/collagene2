//! Target-free reciprocal extension retains literal products and exact ancestry.
use collagene_core::{
    chemistry::{FivePrime, Strand, ThreePrime},
    cloning::overlap::extension::{
        self as overlap_extension,
        pcr::{PrimerPair, Request as Preparation},
    },
    documents::workspace::{self, Workspace},
    material::{
        derivation as native,
        lineage::{EventBinding, InputEventBinding},
        run::{self as run_material, MaterialOrigin},
    },
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    sequence::analysis::Primer,
};
use collagene_model::{Feature, Location, Molecule, Topology};
const ORIGINAL: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTAGATCGTCA";
const EDITED: &str = "AGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCA";
const AMPLICON: &str = "ATGCGAATTCAGTCCGATGCTACGTATTTTCAGACCTGTAGATCGTCAGAATTCGCAT";
fn pair(label: &str, f: &str, ft: &str, r: &str, rt: &str) -> PrimerPair {
    PrimerPair {
        forward: Primer::new(format!("{label} F"), f, ft).unwrap(),
        reverse: Primer::new(format!("{label} R"), r, rt).unwrap(),
        template_strand: Strand::Top,
    }
}
fn preparation() -> Preparation {
    Preparation {
        left: pair("New left", "AGTCCGAT", "", "TACGTAGC", "GGTCTGAAAA"),
        right: pair("New right", "TTCAGACC", "GCTACGTATT", "TGACGATC", ""),
        outer: pair(
            "New outer",
            "AGTCCGAT",
            "ATGCGAATTC",
            "TGACGATC",
            "ATGCGAATTC",
        ),
        extension: overlap_extension::Options {
            overlap: 18,
            min_paired_bases: 8,
            ..Default::default()
        },
    }
}
fn pin(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn original() -> Workspace {
    let mut source = Molecule::dna("Immutable original40", ORIGINAL, Topology::Linear).unwrap();
    source.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(0, 4),
        qualifiers: vec![("label".into(), Some("left source mark".into()))],
        needs_review: false,
    });
    Workspace::from_documents("Target-free source", vec![source], 0).unwrap()
}
fn request(
    w: &Workspace,
    source: EventBinding,
    origin: native::SourceOrigin,
    preparation: Preparation,
) -> native::Request {
    let preview =
        native::overlap_extension::preview_overlap_extension(w, &source, &origin, &preparation)
            .unwrap();
    native::Request {
        sources: vec![collagene_core::material::derivation::BoundSource {
            role: collagene_core::material::derivation::SourceRole::Template,
            event: source,
            origin: origin,
        }],
        derivation: native::Derivation::OverlapExtensionPcr {
            preparation,
            selection: preview.selection().unwrap(),
        },
        name: "Standalone outer-PCR material".into(),
        save_primers: true,
    }
}
#[test]
fn literal_edited_tailed_outer_product_and_all_six_new_oligos_survive_restart_and_digest() {
    let w = original();
    let source = pin(&w, 0);
    let before = workspace::save_workspace(&w).unwrap();
    let prep = preparation();
    let preview = native::overlap_extension::preview_overlap_extension(
        &w,
        &source,
        &native::SourceOrigin::OriginalReference,
        &prep,
    )
    .unwrap();
    preview
        .validate_binding(
            &native::overlap_extension::overlap_extension_preview_revision(
                &w,
                &source,
                &native::SourceOrigin::OriginalReference,
                &prep,
            )
            .unwrap(),
        )
        .unwrap();
    assert!(preview.reusable());
    assert!(!preview.incomplete);
    assert_eq!(
        preview.extension.as_ref().unwrap().fill.duplex.top.as_str(),
        EDITED
    );
    assert_eq!(
        preview.product.as_ref().unwrap().sequence.as_str(),
        AMPLICON
    );
    assert_eq!(
        preview
            .pcr_stages
            .iter()
            .map(|p| p.product_length.unwrap())
            .collect::<Vec<_>>(),
        vec![26, 30, 58]
    );
    assert!(
        preview
            .product
            .as_ref()
            .unwrap()
            .features
            .iter()
            .any(|f| f.label() == "left source mark" && f.location == Location::range(10, 14))
    );
    let req = request(
        &w,
        source.clone(),
        native::SourceOrigin::OriginalReference,
        prep,
    );
    let installed = native::install(&w, &req).unwrap();
    assert_eq!(workspace::save_workspace(&w).unwrap(), before);
    let reopened = workspace::load_workspace(&installed.archive).unwrap();
    installed.receipt.validate_archive(&reopened).unwrap();
    let fresh = native::reconstruct(&reopened, &installed.receipt.product).unwrap();
    assert_eq!(fresh.reference().name, "Standalone outer-PCR material");
    assert_eq!(fresh.reference().sequence.as_str(), AMPLICON);
    assert_eq!(fresh.tracked().duplex.top.as_str(), AMPLICON);
    assert!(fresh.tracked().duplex.nicks.is_empty());
    for terminus in [fresh.tracked().duplex.left, fresh.tracked().duplex.right] {
        assert_eq!(terminus.five_prime, FivePrime::Hydroxyl);
        assert_eq!(terminus.three_prime, ThreePrime::Hydroxyl);
    }
    let asset = &reopened.attachments[installed.receipt.primer_attachment.as_ref().unwrap()];
    let table: collagene_core::formats::oligos::OligoTable =
        serde_json::from_slice(&asset.bytes).unwrap();
    assert_eq!(table.records.len(), 6);
    assert!(table.records.iter().all(|r| r.name.starts_with("New ")));
    let capsule = native::capsule(&reopened, &installed.receipt.product).unwrap();
    let record = native::load_record(&capsule.record_bytes).unwrap();
    assert_eq!(
        native::primer_table(&record)
            .unwrap()
            .unwrap()
            .records
            .len(),
        6
    );
    let catalog = collagene_core::cloning::restriction::catalog::builtin_catalog().unwrap();
    let selected =
        collagene_core::cloning::restriction::catalog::select(&catalog, &["ecori".into()]).unwrap();
    let digest = native::digest::inspect_digest(
        "outer-product",
        fresh.reference(),
        &native::SourceOrigin::Saved(MaterialOrigin::Derived(capsule)),
        &selected,
        &Default::default(),
    )
    .unwrap();
    assert!(digest.inspection.products_determined);
    assert_eq!(
        digest
            .inspection
            .fragments
            .iter()
            .map(|f| f.duplex.top_length)
            .collect::<Vec<_>>(),
        vec![5, 44, 9]
    );
}
#[test]
fn competing_pcrs_withhold_material_and_stale_or_altered_previews_do_not_authorize_install() {
    let w = original();
    let source = pin(&w, 0);
    let prep = preparation();
    let origin = native::SourceOrigin::OriginalReference;
    let preview =
        native::overlap_extension::preview_overlap_extension(&w, &source, &origin, &prep).unwrap();
    let mut altered = preview.clone();
    altered.preparation.outer.forward.tail.push('A');
    assert!(altered.validate_binding(&preview.context_revision).is_err());
    let mut altered = preview.clone();
    altered.product.as_mut().unwrap().name.push('!');
    assert!(altered.validate_binding(&preview.context_revision).is_err());
    let mut req = request(&w, source.clone(), origin.clone(), prep.clone());
    if let native::Derivation::OverlapExtensionPcr { selection, .. } = &mut req.derivation {
        selection.expected_product_revision = "a".repeat(64)
    }
    assert!(native::install(&w, &req).is_err());
    let repeated = Workspace::from_documents(
        "Repeated sites",
        vec![
            Molecule::dna(
                "Repeated",
                &format!("{ORIGINAL}{ORIGINAL}"),
                Topology::Linear,
            )
            .unwrap(),
        ],
        0,
    )
    .unwrap();
    let preview = native::overlap_extension::preview_overlap_extension(
        &repeated,
        &pin(&repeated, 0),
        &origin,
        &prep,
    )
    .unwrap();
    assert!(preview.incomplete);
    assert!(preview.product.is_none());
    assert!(preview.selection().is_err());
    assert!(preview.pcr_stages[0].geometries.len() > 1);
    let mut prep = prep;
    prep.outer = pair("Competing outer", "T", "", "A", "");
    let preview =
        native::overlap_extension::preview_overlap_extension(&w, &source, &origin, &prep).unwrap();
    assert!(preview.incomplete);
    assert!(preview.product.is_none());
    assert!(preview.selection().is_err());
}
fn saved_run(incomplete: bool) -> (Workspace, EventBinding, MaterialOrigin) {
    let source = Molecule::dna("Origin source", ORIGINAL, Topology::Linear).unwrap();
    let mut steps = vec![];
    let target = if incomplete {
        let mut p = preparation();
        p.outer.forward.tail.clear();
        p.outer.reverse.tail.clear();
        p.extension.max_comparisons = 18;
        steps = vec![
            Operation::Pcr {
                template: Source::Input("source".into()),
                forward: p.left.forward,
                reverse: p.left.reverse,
                template_strand: Strand::Top,
            },
            Operation::Pcr {
                template: Source::Input("source".into()),
                forward: p.right.forward,
                reverse: p.right.reverse,
                template_strand: Strand::Top,
            },
            Operation::OverlapExtension {
                left: collagene_core::planning::plan::OrientedSource {
                    source: Source::Step(0),
                    reverse: false,
                },
                right: collagene_core::planning::plan::OrientedSource {
                    source: Source::Step(1),
                    reverse: false,
                },
                options: p.extension,
            },
            Operation::Pcr {
                template: Source::Step(2),
                forward: p.outer.forward,
                reverse: p.outer.reverse,
                template_strand: Strand::Top,
            },
        ];
        Molecule::dna("Independent source-run target", EDITED, Topology::Linear).unwrap()
    } else {
        steps.push(Operation::Pcr {
            template: Source::Input("source".into()),
            forward: Primer::new("Old F", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("Old R", "TGACGATC", "").unwrap(),
            template_strand: Strand::Top,
        });
        Molecule::dna("Independent source-run target", ORIGINAL, Topology::Linear).unwrap()
    };
    let plan = DesignPlan {
        name: "Original source amplification".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "source".into(),
            revision: source.revision(),
            molecule: source.clone(),
        }],
        steps,
    };
    let w = Workspace::from_documents("Exact originating run", vec![source, target], 1).unwrap();
    let commit = collagene_core::material::outputs::finalize(
        &w,
        &plan,
        pin(&w, 1),
        vec![InputEventBinding {
            input_id: "source".into(),
            event: pin(&w, 0),
        }],
        &collagene_core::material::outputs::Options {
            name: "Source run product".into(),
            show_prepared: false,
            save_primers: true,
            allow_incomplete: incomplete,
        },
    )
    .unwrap();
    let event = commit.receipt.construct;
    let origin = collagene_core::material::capsule::capsule(
        &commit.workspace,
        &run_material::Selection {
            run_id: commit.receipt.run.run_id,
            step_index: plan.steps.len() - 1,
            event: event.clone(),
        },
    )
    .unwrap();
    (commit.workspace, event, MaterialOrigin::Run(origin))
}
#[test]
fn saved_run_and_derived_origins_keep_ownership_and_do_not_reexport_historical_primers() {
    let (w, source, origin) = saved_run(false);
    let req = request(
        &w,
        source.clone(),
        native::SourceOrigin::Saved(origin.clone()),
        preparation(),
    );
    let foreign = original();
    assert!(
        native::overlap_extension::preview_overlap_extension(
            &foreign,
            &source,
            &native::SourceOrigin::Saved(origin.clone()),
            &preparation()
        )
        .is_err()
    );
    assert!(
        native::overlap_extension::preview_overlap_extension(
            &w,
            &source,
            &native::SourceOrigin::OriginalReference,
            &preparation()
        )
        .is_err()
    );
    let mut alias = source.clone();
    alias.event_id = pin(&w, 0).event_id;
    assert!(
        native::overlap_extension::preview_overlap_extension(
            &w,
            &alias,
            &native::SourceOrigin::Saved(origin),
            &preparation()
        )
        .is_err()
    );
    let first = native::prepare_install(&w, &req).unwrap();
    let new_table: collagene_core::formats::oligos::OligoTable = serde_json::from_slice(
        &first.workspace.attachments[first.receipt.primer_attachment.as_ref().unwrap()].bytes,
    )
    .unwrap();
    assert_eq!(new_table.records.len(), 6);
    assert!(new_table.records.iter().all(|r| r.name.starts_with("New ")));
    let source = first.receipt.product;
    let origin = MaterialOrigin::Derived(native::capsule(&first.workspace, &source).unwrap());
    let req = request(
        &first.workspace,
        source.clone(),
        native::SourceOrigin::Saved(origin),
        preparation(),
    );
    let second = native::install(&first.workspace, &req).unwrap();
    let reopened = workspace::load_workspace(&second.archive).unwrap();
    assert_eq!(
        native::reconstruct(&reopened, &second.receipt.product)
            .unwrap()
            .reference()
            .sequence
            .as_str(),
        AMPLICON
    );
    let request_source = first.workspace.revisions[&source.revision].clone();
    assert!(
        native::overlap_extension::preview_overlap_extension(
            &original(),
            &source,
            &native::SourceOrigin::OriginalReference,
            &preparation()
        )
        .is_err()
    );
    assert_eq!(request_source.sequence.as_str(), AMPLICON);
}
#[test]
fn inherited_origin_and_terminal_search_limits_remain_conditional_after_selection_and_restart() {
    let (w, source, origin) = saved_run(true);
    let req = request(
        &w,
        source.clone(),
        native::SourceOrigin::Saved(origin),
        preparation(),
    );
    let preview = native::overlap_extension::preview_overlap_extension(
        &w,
        &source,
        &req.sources[0].origin,
        &preparation(),
    )
    .unwrap();
    assert!(preview.origin_incomplete);
    assert!(preview.incomplete);
    assert!(preview.reusable());
    let mut changed = preview.clone();
    changed.incomplete = false;
    assert!(changed.validate_binding(&preview.context_revision).is_err());
    let installed = native::install(&w, &req).unwrap();
    let reopened = workspace::load_workspace(&installed.archive).unwrap();
    assert!(installed.receipt.incomplete);
    assert!(
        native::reconstruct(&reopened, &installed.receipt.product)
            .unwrap()
            .incomplete()
    );
    let w = original();
    let mut prep = preparation();
    prep.extension.max_comparisons = 18;
    let req = request(
        &w,
        pin(&w, 0),
        native::SourceOrigin::OriginalReference,
        prep,
    );
    let installed = native::install(&w, &req).unwrap();
    assert!(installed.receipt.incomplete);
    assert_eq!(
        native::reconstruct(
            &workspace::load_workspace(&installed.archive).unwrap(),
            &installed.receipt.product
        )
        .unwrap()
        .reference()
        .sequence
        .as_str(),
        AMPLICON
    );
}
