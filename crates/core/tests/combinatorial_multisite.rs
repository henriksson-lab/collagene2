//! Whole original-parent reference networks in the existing selected-row compiler.
use collagene_core::{
    Molecule, Topology, assembly_recipe, chemistry, cloning_outputs, combinatorial_batch as batch,
    combinatorial_design::{
        self as design, Bin, Choice, ChoiceSelection, ConstructOverrides, Design, Method,
        MultiSiteSelection, PreparationOverride, SelectedConstruct, TargetDeclaration,
    },
    combinatorial_export as export, combinatorial_release as release,
    multisite_recombination::{self as multi, multisite},
    recombination::reference::{Family, Reaction, SiteVariant, SourceRecord},
    run_lineage::EventBinding,
    verification::{self, Operation, Source, Status},
    workspace::{self, Workspace},
};
const MAIN: &str = concat!(
    "GTACAAT", "TCC", "CCCCA", "TGGA", "GTACAAA", "TCC", "GGGGT", "ACGG", "GTACAAG", "GCT",
    "AAAAG", "ACGG", "GTACAAC", "GCT", "GGGGA", "TGGA"
);
fn core(id: u8) -> &'static str {
    match id {
        1 => "GTACAAA",
        2 => "GTACAAG",
        3 => "GTACAAC",
        4 => "GTACAAT",
        _ => panic!("fixture specificity"),
    }
}
fn site(family: Family, id: u8) -> String {
    match family {
        Family::L => format!("ACGG{}TCC", core(id)),
        Family::R => format!("TGGA{}GCT", core(id)),
        _ => panic!("fixture family"),
    }
}
fn event(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn fixture(existing: bool) -> (Workspace, Design, SelectedConstruct) {
    let mut profile=multisite::Profile{id:"synthetic-library-four-parents".into(),edition:"literal-source-fixture-v1".into(),sources:vec![],variants:vec![],limitations:vec!["Explicit synthetic reference fixture; no empirical or manufacturer-specific applicability claim".into()]};
    for id in 1..=4 {
        for (family, prefix) in [(Family::L, "L"), (Family::R, "R")] {
            let key = format!("{prefix}{id}");
            let sequence = site(family, id);
            profile.sources.push(SourceRecord {
                id: key.clone(),
                url: "https://example.org/synthetic-fixture".into(),
                raw_sha256: "1".repeat(64),
                vector_name: "Synthetic declared reference context".into(),
                vector_length: sequence.len(),
                start_1based: 1,
                end_1based: sequence.len(),
                reverse_extraction: false,
                map_url: "https://example.org/synthetic-fixture-map".into(),
                map_sha256: "2".repeat(64),
            });
            profile.variants.push(SiteVariant {
                id: key.clone(),
                family,
                subtype: id,
                sequence,
                overlap_start: 4,
                source_id: key,
            });
        }
    }
    let strings = [
        format!("{}CCCCA{}AAAAC", site(Family::L, 4), site(Family::R, 1)),
        format!("{}GGGGT{}TTTTC", site(Family::L, 1), site(Family::L, 2)),
        format!("{}AAAAG{}CCCCT", site(Family::R, 2), site(Family::L, 3)),
        format!("{}TTTTA{}GGGGA", site(Family::R, 4), site(Family::R, 3)),
    ];
    let mut docs = strings
        .iter()
        .enumerate()
        .map(|(i, s)| Molecule::dna(format!("Original parent {i}"), s, Topology::Circular).unwrap())
        .collect::<Vec<_>>();
    docs.push(Molecule::dna("Separate immutable target", MAIN, Topology::Circular).unwrap());
    docs.push(docs[0].clone());
    docs.push(docs[4].clone());
    docs.push(
        Molecule::dna(
            "Wrong independent target",
            &"A".repeat(76),
            Topology::Circular,
        )
        .unwrap(),
    );
    let w = Workspace::from_documents("Whole original multi-site parent choices", docs, 4).unwrap();
    let regions = [(4, 23), (4, 23), (4, 23), (23, 42)];
    let bins = regions
        .iter()
        .enumerate()
        .map(|(i, (start, end))| Bin {
            id: format!("parent-{i}"),
            label: format!("Parent {i}"),
            choices: vec![Choice {
                generation: None,
                id: format!("parent-{i}-choice"),
                label: format!("Original parent {i} with explicit intended contribution"),
                source: event(&w, i),
                start: *start,
                end: *end,
                reverse: false,
                preparation: None,
                treatment: None,
            }],
        })
        .collect::<Vec<_>>();
    let row = SelectedConstruct {
        id: "wanted-reference".into(),
        name: "User chosen clone name".into(),
        choices: bins
            .iter()
            .map(|bin| ChoiceSelection {
                bin_id: bin.id.clone(),
                choice_id: bin.choices[0].id.clone(),
            })
            .collect(),
        target: if existing {
            TargetDeclaration::Existing(event(&w, 4))
        } else {
            TargetDeclaration::OriginalRegions {
                junction_edits: vec![],
            }
        },
        overrides: None,
    };
    let d = Design {
        id: "whole-parent-library".into(),
        name: "Connected multi-site library".into(),
        circular: true,
        bins,
        method: Method::MultiSiteReference {
            request: multi::Request {
                reaction: multisite::Request {
                    profile,
                    reaction: Reaction::LR,
                    specificities: vec![4, 1, 2, 3],
                    limits: Default::default(),
                },
                annotation_limits: Default::default(),
            },
            selection: None,
        },
    };
    (w, d, row)
}
fn pin(w: &Workspace, d: &mut Design, row: &SelectedConstruct) -> design::MultiSiteInspection {
    let inspected = design::inspect_multisite(d, row, w).unwrap();
    // The main product contains this original interval forward; do not choose a
    // component by searching the independently supplied target sequence.
    let component = inspected.batch.reference.candidates[0]
        .products
        .iter()
        .position(|product| {
            product.segments.iter().any(|map| {
                map.input == 0 && map.source_start == 4 && map.source_end == 23 && !map.reverse
            })
        })
        .unwrap();
    let Method::MultiSiteReference { selection, .. } = &mut d.method else {
        panic!("fixture method")
    };
    *selection = Some(MultiSiteSelection {
        selection: multi::Selection {
            batch_revision: inspected.batch.batch_revision.clone(),
            candidate_id: 0,
            component_id: component,
        },
        source_revision: inspected.source_revision.clone(),
    });
    inspected
}
#[test]
fn original_contribution_intent_is_frozen_before_whole_parent_discovery_and_explicit_selection() {
    let (w, mut d, row) = fixture(false);
    assert_eq!(postcard::to_allocvec(&d.method).unwrap()[0], 4);
    let resolved = design::resolve_multisite(&d, &row, &w).unwrap();
    assert_eq!(resolved.target.sequence.as_str(), MAIN);
    assert!(resolved.target_binding.is_none());
    assert!(
        resolved
            .inputs
            .iter()
            .all(|input| input.molecule.sequence.len() == 38)
    );
    assert!(design::compile(&d, &row, &w).is_err());
    let before_context = design::inspect_multisite_revision(&d, &row, &w).unwrap();
    let inspected = pin(&w, &mut d, &row);
    assert_eq!(inspected.context_revision, before_context);
    assert_eq!(inspected.batch.documents.len(), 2);
    assert_eq!(
        inspected.batch.reference.candidates[0]
            .products
            .iter()
            .map(|p| p.molecule.sequence.len())
            .sum::<usize>(),
        152
    );
    let after = design::resolve_multisite(&d, &row, &w).unwrap();
    let design::Recipe::MultiSiteReference {
        source_revision, ..
    } = &after.recipe
    else {
        panic!("reference recipe")
    };
    assert_eq!(source_revision, &inspected.source_revision);
    let compiled = design::compile(&d, &row, &w).unwrap();
    assert_eq!(compiled.plan.name, row.name);
    assert!(
        matches!(&compiled.plan.steps[..],[Operation::MultiSiteRecombineReference{sources,..}] if sources.len()==4&&sources.iter().all(|s|matches!(s.source,Source::Input(_))))
    );
    let evaluated = design::evaluate(&d, &row, &w).unwrap();
    assert_eq!(evaluated.verification().status, Status::Incomplete);
    assert!(evaluated.preparation_complete());
    assert_eq!(evaluated.verification().products.len(), 1);
    assert_eq!(evaluated.verification().reference_step_reports.len(), 1);
    assert!(chemistry::Duplex::from_molecule(&evaluated.verification().products[0]).is_err());
    let oligos = collagene_core::plan_oligos::build(
        evaluated.plan(),
        &evaluated.resolved.target,
        &Default::default(),
    )
    .unwrap();
    assert!(oligos.uses.is_empty());
    assert!(oligos.oligos.records.is_empty());
    let reinspection = design::inspect_multisite(&d, &row, &w).unwrap();
    assert_eq!(reinspection.source_revision, inspected.source_revision);
    assert_eq!(
        reinspection.batch.batch_revision,
        inspected.batch.batch_revision
    );
    assert_ne!(reinspection.context_revision, inspected.context_revision);
}
#[test]
fn exact_event_or_target_changes_invalidate_pins_even_when_literal_dna_matches() {
    let (w, mut d, mut row) = fixture(true);
    let inspected = pin(&w, &mut d, &row);
    assert_eq!(design::compile(&d, &row, &w).unwrap().plan.name, row.name);
    let mut same_dna_source = d.clone();
    same_dna_source.bins[0].choices[0].source = event(&w, 5);
    assert_eq!(event(&w, 0).revision, event(&w, 5).revision);
    assert!(design::compile(&same_dna_source, &row, &w).is_err());
    let moved = design::inspect_multisite(&same_dna_source, &row, &w).unwrap();
    assert_eq!(moved.batch.batch_revision, inspected.batch.batch_revision);
    assert_ne!(moved.source_revision, inspected.source_revision);
    row.target = TargetDeclaration::Existing(event(&w, 6));
    assert_eq!(event(&w, 4).revision, event(&w, 6).revision);
    assert!(design::compile(&d, &row, &w).is_err());
    let moved_target = design::inspect_multisite(&d, &row, &w).unwrap();
    assert_eq!(
        moved_target.batch.batch_revision,
        inspected.batch.batch_revision
    );
    assert_ne!(moved_target.source_revision, inspected.source_revision);
    row.target = TargetDeclaration::Existing(event(&w, 7));
    assert!(design::compile(&d, &row, &w).is_err());
    let wrong = pin(&w, &mut d, &row);
    assert_eq!(wrong.batch.batch_revision, inspected.batch.batch_revision);
    let evaluated = design::evaluate(&d, &row, &w).unwrap();
    assert_eq!(evaluated.verification().status, Status::Invalid);
    assert_eq!(evaluated.verification().products.len(), 1);
    assert_eq!(evaluated.verification().reference_step_reports.len(), 1);
}
#[test]
fn preparation_treatment_and_topology_are_explicit_and_per_construct_reference_overrides_work() {
    let (w, mut d, mut row) = fixture(true);
    d.bins[0].choices[0].preparation = Some(assembly_recipe::Preparation::WholeLinear);
    d.bins[0].choices[0].treatment = Some(chemistry::EndTreatment::Phosphorylate);
    assert!(design::inspect_multisite(&d, &row, &w).is_err());
    row.overrides = Some(ConstructOverrides {
        preparations: vec![PreparationOverride {
            generation: None,
            bin_id: d.bins[0].id.clone(),
            preparation: None,
            treatment: None,
        }],
        method: None,
        circular: None,
    });
    pin(&w, &mut d, &row);
    assert_eq!(
        design::evaluate(&d, &row, &w)
            .unwrap()
            .verification()
            .status,
        Status::Incomplete
    );
    let multisite = d.method.clone();
    d.method = Method::Ligation {
        allow_nicked: false,
    };
    row.overrides.as_mut().unwrap().method = Some(multisite);
    assert_eq!(
        design::evaluate(&d, &row, &w)
            .unwrap()
            .verification()
            .status,
        Status::Incomplete
    );
    row.overrides.as_mut().unwrap().circular = Some(false);
    assert!(design::inspect_multisite(&d, &row, &w).is_err());
}
#[test]
fn existing_batch_release_and_export_preserve_incomplete_reference_status_and_all_graphs() {
    let (w, mut d, row) = fixture(false);
    pin(&w, &mut d, &row);
    let mut session = batch::Session::begin(
        w,
        batch::Request {
            design: Some(d),
            rows: vec![batch::RowRequest::Ready(row.clone())],
        },
    )
    .unwrap();
    let context = session.summary().context_revision;
    let progress = session.advance(&context, 0).unwrap();
    let evidence = progress.row.unwrap();
    assert_eq!(evidence.status, Some(Status::Incomplete));
    assert!(evidence.preparation_complete);
    let captured = session.into_capture(&context).unwrap();
    let mut options = cloning_outputs::Options {
        name: "Released reference clone".into(),
        ..Default::default()
    };
    let request = release::Request {
        batch_hash: captured.attachment_hash.clone(),
        rows: vec![release::Selection {
            construct_id: row.id.clone(),
            options: options.clone(),
        }],
    };
    let mut refused = release::ReleaseSession::begin(captured.workspace.clone(), request).unwrap();
    let failed_context = refused.summary().context_revision;
    assert!(refused.advance(&failed_context, 0).is_err());
    assert!(refused.into_commit(&failed_context).is_err());
    options.allow_incomplete = true;
    let mut released = release::ReleaseSession::begin(
        captured.workspace,
        release::Request {
            batch_hash: captured.attachment_hash,
            rows: vec![release::Selection {
                construct_id: row.id.clone(),
                options,
            }],
        },
    )
    .unwrap();
    let release_context = released.summary().context_revision;
    released.advance(&release_context, 0).unwrap();
    let inspection = released.inspect(0).unwrap();
    assert!(chemistry::has_unmapped_physical_evidence(
        &inspection.product
    ));
    assert_eq!(
        inspection.row.output.run.recorded_status,
        Status::Incomplete
    );
    let committed = released.into_commit(&release_context).unwrap();
    let manifest = release::load(&committed.workspace, &committed.manifest_hash).unwrap();
    assert_eq!(manifest.rows.len(), 1);
    assert!(manifest.rows[0].options.allow_incomplete);
    let graphs = collagene_core::reference_run_evidence::load(
        &committed.workspace,
        &manifest.rows[0].output.run.run_id,
    )
    .unwrap()
    .unwrap();
    assert_eq!(graphs.steps.len(), 1);
    let verification::ReferenceStepEvidence::MultiSiteRecombination { batch, .. } =
        &graphs.steps[0].evidence;
    assert_eq!(batch.documents.len(), 2);
    assert_eq!(
        batch
            .documents
            .iter()
            .map(|d| d.molecule.sequence.len())
            .sum::<usize>(),
        152
    );
    let reopened =
        workspace::load_workspace(&workspace::save_workspace(&committed.workspace).unwrap())
            .unwrap();
    reopened.validate().unwrap();
    let archive = export::export(
        &reopened,
        &committed.manifest_hash,
        &export::Options {
            genbank: true,
            fasta: true,
            primer_csv: false,
            primer_tsv: false,
            protocol_html: true,
            protocol_pdf: true,
            construct_ids: None,
        },
    )
    .unwrap();
    assert!(archive.bytes.starts_with(b"PK"));
    assert!(archive.files.iter().any(|f| f.path.ends_with(".html")));
    assert!(archive.files.iter().any(|f| f.path.ends_with(".pdf")));
    assert!(archive.files.iter().any(|f| f.path.ends_with(".gb")));
}
