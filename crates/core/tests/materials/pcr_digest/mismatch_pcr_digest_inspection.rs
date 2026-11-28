use collagene_core::{
    chemistry::{FivePrime, ThreePrime},
    cloning::{pcr::mismatch_digest as m, restriction::catalog as catalog_digest},
    sequence::analysis::Primer,
};
use collagene_model::{Feature, Location, Molecule, Topology};
fn fixture() -> (Molecule, m::Request) {
    let mut source =
        Molecule::dna("Original template", "ACGACCCCGGGGTAAC", Topology::Linear).unwrap();
    source.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::range(8, 12),
        qualifiers: vec![("label".into(), Some("retained domain".into()))],
        needs_review: false,
    });
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let id = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRI")
        .unwrap()
        .id
        .clone();
    (
        source,
        m::Request {
            base: collagene_core::cloning::pcr::digest::Request {
                source_role: "original".into(),
                forward: Primer::new("Mutation F", "ACTACCCC", "CCGAATTC").unwrap(),
                reverse: Primer::new("Mutation R", "GTTACCCC", "GGGAATTC").unwrap(),
                selection: catalog_digest::select(&catalog, &[id]).unwrap(),
                options: Default::default(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            pcr_options: Default::default(),
        },
    )
}
#[test]
fn literal_mismatch_tail_digest_and_annotation_positions_retain_actual_cut_chemistry() {
    let (source, request) = fixture();
    let before = serde_json::to_value(&source).unwrap();
    let result = m::inspect(&source, None, &request, None).unwrap();
    result
        .validate_binding(&m::inspection_revision(&source, None, &request, None).unwrap())
        .unwrap();
    assert_eq!(result.pcr.report.products.len(), 1);
    assert!(result.pcr.report.enumeration_complete);
    let product = result.amplicon.as_ref().unwrap();
    assert_eq!(
        product.sequence.as_str(),
        "CCGAATTCACTACCCCGGGGTAACGAATTCCC"
    );
    assert_eq!(
        result.pcr.report.products[0].substitutions[0].source_position,
        2
    );
    assert_eq!(
        result.pcr.report.products[0].substitutions[0].product_position,
        10
    );
    assert!(product.features.iter().any(|feature| {
        feature.label() == "retained domain"
            && feature
                .location
                .spans()
                .iter()
                .any(|span| span.0 == 16 && span.1 == 20)
    }));
    let report = result.inspection.as_ref().unwrap();
    assert!(report.scan_complete && report.products_determined);
    assert_eq!(report.fragments.len(), 3);
    let ends = report
        .fragments
        .iter()
        .map(|f| f.duplex.ends.as_ref().unwrap())
        .collect::<Vec<_>>();
    assert_eq!(ends[0].0.chemistry.five_prime, FivePrime::Hydroxyl);
    assert_eq!(ends[0].1.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(ends[1].0.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(ends[1].1.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(ends[2].0.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(ends[2].1.chemistry.five_prime, FivePrime::Hydroxyl);
    for (left, right) in ends {
        assert_eq!(left.chemistry.three_prime, ThreePrime::Hydroxyl);
        assert_eq!(right.chemistry.three_prime, ThreePrime::Hydroxyl);
    }
    assert!(!result.incomplete);
    assert_eq!(serde_json::to_value(&source).unwrap(), before);
}
#[test]
fn partial_singleton_is_explicitly_conditional_but_multiple_zero_and_overlap_withhold_choices() {
    let (_, mut request) = fixture();
    let source = Molecule::dna(
        "Repeated original",
        "ACGACCCCGGGGTAACACGACCCCGGGGTAAC",
        Topology::Linear,
    )
    .unwrap();
    let multiple = m::inspect(&source, None, &request, None).unwrap();
    assert!(multiple.pcr.report.products.len() > 1);
    assert!(multiple.amplicon.is_none() && multiple.inspection.is_none());
    request.pcr_options.max_products = 1;
    let partial = m::inspect(&source, None, &request, None).unwrap();
    assert_eq!(partial.pcr.report.products.len(), 1);
    assert!(!partial.pcr.report.enumeration_complete);
    assert!(partial.incomplete);
    assert!(partial.amplicon.is_some() && partial.inspection.is_some());
    let empty = Molecule::dna("No placements", "TTTTTTTTTTTTTTTT", Topology::Linear).unwrap();
    let zero = m::inspect(&empty, None, &request, None).unwrap();
    assert!(zero.pcr.report.products.is_empty());
    assert!(zero.amplicon.is_none() && zero.inspection.is_none());
    let source = Molecule::dna("Overlapping domains", "ACGTAC", Topology::Linear).unwrap();
    request.base.forward = Primer::new("F", "ACGT", "").unwrap();
    request.base.reverse = Primer::new("R", "GTAC", "").unwrap();
    request.pcr_options.binding.max_mismatches = 0;
    request.pcr_options.max_products = 1000;
    let overlap = m::inspect(&source, None, &request, None).unwrap();
    assert!(
        overlap
            .pcr
            .report
            .products
            .iter()
            .any(|p| p.overlapping_domains)
    );
    assert!(overlap.amplicon.is_none() && overlap.inspection.is_none());
}
#[test]
fn readonly_preview_binding_rejects_source_policy_catalog_document_and_completion_tampering() {
    let (source, request) = fixture();
    let original = m::inspect(&source, None, &request, None).unwrap();
    let mut changed = request.clone();
    changed.pcr_options.binding.max_mismatches = 0;
    assert!(
        original
            .validate_binding(&m::inspection_revision(&source, None, &changed, None).unwrap())
            .is_err()
    );
    let mut changed = original.clone();
    changed.pcr.report.bindings[0].annealing = "TTTT".into();
    assert!(
        changed
            .validate_binding(&original.context_revision)
            .is_err()
    );
    let mut changed = original.clone();
    changed.pcr.documents[0]
        .metadata
        .insert("pcr:template_revision".into(), "0".repeat(64));
    assert!(
        changed
            .validate_binding(&original.context_revision)
            .is_err()
    );
    let mut changed = original.clone();
    changed.request.base.selection.catalog_content_hash = "0".repeat(64);
    assert!(
        changed
            .validate_binding(&original.context_revision)
            .is_err()
    );
    let mut changed = original.clone();
    changed.incomplete = true;
    assert!(
        changed
            .validate_binding(&original.context_revision)
            .is_err()
    );
    let mut marked = source;
    marked.metadata.insert("duplex_model".into(), "{}".into());
    assert!(m::source_preflight(&marked, None, &request, None).is_err());
}
#[test]
fn saved_origin_is_freshly_replayed_and_forged_origin_or_stale_source_is_not_authority() {
    use collagene_core::{
        documents::workspace::Workspace,
        material::{
            lineage::{EventBinding, InputEventBinding},
            outputs as cloning_outputs, run as run_material,
        },
        planning::plan::{DesignPlan, Operation, PlanInput, Source},
    };
    let (source, request) = fixture();
    let target = Molecule::dna(
        "Independent original PCR",
        "ACGACCCCGGGGTAAC",
        Topology::Linear,
    )
    .unwrap();
    let workspace =
        Workspace::from_documents("Origin", vec![source.clone(), target.clone()], 0).unwrap();
    let event = |i: usize| {
        let d = &workspace.documents[i];
        EventBinding {
            document_id: d.id.clone(),
            event_id: d.head_event.clone().unwrap(),
            revision: d.head.clone(),
        }
    };
    let plan = DesignPlan {
        name: "Original exact PCR".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "original".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("original".into()),
            forward: Primer::new("Origin F", "ACGACCCC", "").unwrap(),
            reverse: Primer::new("Origin R", "GTTACCCC", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        }],
    };
    let committed = cloning_outputs::finalize(
        &workspace,
        &plan,
        event(1),
        vec![InputEventBinding {
            input_id: "original".into(),
            event: event(0),
        }],
        &cloning_outputs::Options {
            name: "Saved exact product".into(),
            ..Default::default()
        },
    )
    .unwrap();
    let source = committed
        .workspace
        .head(&committed.receipt.construct.document_id)
        .unwrap();
    let origin = collagene_core::material::capsule::capsule(
        &committed.workspace,
        &run_material::Selection {
            run_id: committed.receipt.run.run_id,
            step_index: 0,
            event: committed.receipt.construct,
        },
    )
    .unwrap();
    let result = m::inspect(
        source,
        Some(&run_material::MaterialOrigin::Run(origin.clone())),
        &request,
        Some(&committed.workspace),
    )
    .unwrap();
    assert!(result.origin_revision.is_some() && result.origin_incomplete == Some(false));
    assert_eq!(
        result.amplicon.unwrap().sequence.as_str(),
        "CCGAATTCACTACCCCGGGGTAACGAATTCCC"
    );
    let mut forged = origin.clone();
    forged.products[0].revision = "0".repeat(64);
    assert!(
        m::inspect(
            source,
            Some(&run_material::MaterialOrigin::Run(forged.clone())),
            &request,
            Some(&committed.workspace)
        )
        .is_err()
    );
    let mut changed = source.clone();
    changed.name.push_str(" stale");
    assert!(
        m::source_preflight(
            &changed,
            Some(&run_material::MaterialOrigin::Run(origin.clone())),
            &request,
            Some(&committed.workspace)
        )
        .is_err()
    );
    // A genuine target-free digest derivative is retained and reopened before
    // the same explicit mismatch-PCR operation; it is not cast to a Run.
    use collagene_core::{documents::workspace, material::derivation as d};
    let catalog = catalog_digest::builtin_catalog().unwrap();
    let id = catalog
        .edition()
        .enzymes
        .iter()
        .find(|e| e.name == "EcoRV")
        .unwrap()
        .id
        .clone();
    let prepared = d::prepare_install(
        &committed.workspace,
        &d::Request {
            sources: vec![collagene_core::material::derivation::BoundSource {
                role: collagene_core::material::derivation::SourceRole::Template,
                event: origin.selection.event.clone(),
                origin: d::SourceOrigin::Saved(run_material::MaterialOrigin::Run(origin.clone())),
            }],
            derivation: d::Derivation::CatalogDigest {
                selection: catalog_digest::select(&catalog, &[id]).unwrap(),
                options: Default::default(),
                fragment_index: 0,
            },
            name: "Retained native digest template".into(),
            save_primers: false,
        },
    )
    .unwrap();
    let reopened =
        workspace::load_workspace(&workspace::save_workspace(&prepared.workspace).unwrap())
            .unwrap();
    let derived = run_material::MaterialOrigin::Derived(
        d::capsule(&reopened, &prepared.receipt.product).unwrap(),
    );
    let template = &reopened.revisions[&derived.event().revision];
    let before = workspace::save_workspace(&reopened).unwrap();
    let report = m::inspect(template, Some(&derived), &request, Some(&reopened)).unwrap();
    assert_eq!(report.origin_incomplete, Some(false));
    assert_eq!(
        report.amplicon.as_ref().unwrap().sequence.as_str(),
        "CCGAATTCACTACCCCGGGGTAACGAATTCCC"
    );
    let ends = report.inspection.as_ref().unwrap().fragments[1]
        .duplex
        .ends
        .as_ref()
        .unwrap();
    assert_eq!(ends.0.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(ends.1.chemistry.five_prime, FivePrime::Phosphate);
    assert!(report.pcr.report.products[0].certain);
    assert_eq!(
        report.pcr.report.products[0].substitutions[0].source_position,
        2
    );
    assert_eq!(
        report.pcr.report.products[0].substitutions[0].product_position,
        10
    );
    assert_eq!(workspace::save_workspace(&reopened).unwrap(), before);
    assert!(m::inspect(template, Some(&derived), &request, None).is_err());
    assert!(m::inspect(template, Some(&derived), &request, Some(&workspace)).is_err());
    let mut bottom = request.clone();
    bottom.base.template_strand = collagene_core::chemistry::Strand::Bottom;
    assert!(m::inspect(template, Some(&derived), &bottom, Some(&reopened)).is_err());
    let run_material::MaterialOrigin::Derived(mut wrong) = derived.clone() else {
        panic!("Derived")
    };
    wrong.product.event_id.push_str("-alias");
    assert!(
        m::source_preflight(
            template,
            Some(&run_material::MaterialOrigin::Derived(wrong)),
            &request,
            Some(&reopened)
        )
        .is_err()
    );
    let mut changed = request.clone();
    changed.base.forward.tail = "AAGAATTC".into();
    let changed_context =
        m::inspection_revision(template, Some(&derived), &changed, Some(&reopened)).unwrap();
    assert!(report.validate_binding(&changed_context).is_err());
}

#[test]
fn shared_preparation_profile_bounds_ordered_oligos_and_combined_work_before_inspection() {
    let (source, mut request) = fixture();
    // Separate domains/tails pass the generic exact-preview limits, while this
    // ordinary preparation profile bounds each full ordered oligo to200 bases.
    request.base.forward.tail = "A".repeat(193);
    assert!(m::source_preflight(&source, None, &request, None).is_err());
    let (_, mut request) = fixture();
    request.base.forward = Primer::new("Bound F", &"A".repeat(100), "").unwrap();
    request.base.reverse = Primer::new("Bound R", &"C".repeat(100), "").unwrap();
    request.pcr_options.binding.max_work = 200_000_000;
    //100KiB source is small; its declared worst-case mismatch/materialization
    // and catalog work exceeds the shared100M preparation profile nonetheless.
    let source = Molecule::dna("Bounded source", &"A".repeat(100_000), Topology::Linear).unwrap();
    assert!(matches!(
        m::source_preflight(&source, None, &request, None),
        Err(collagene_model::Error::TooLarge)
    ));
}
