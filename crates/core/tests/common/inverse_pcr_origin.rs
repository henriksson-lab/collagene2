use collagene_core::{
    Molecule, Topology,
    analysis::Primer,
    catalog_digest,
    chemistry::EndTreatment,
    cloning_outputs, material_derivation as derive,
    mutagenesis::InversePcrOptions,
    run_lineage::{EventBinding, InputEventBinding},
    run_material::{self, MaterialOrigin},
    verification::{DesignPlan, Operation, OrientedSource, PlanInput, Source},
    workspace::Workspace,
};
pub const STOCK: &str = "AGTCCGATGCTACGTAGGACTTCAGACCTGTA";
pub const INTENDED: &str = "AGTGGAGATGCTACGTAGGACTTCAGACCTGTA";
pub fn options() -> InversePcrOptions {
    InversePcrOptions {
        edit_start: 3,
        edit_end: 5,
        inserted: "GGA".into(),
        forward_length: 8,
        reverse_length: 8,
        reverse_tail_bases: 0,
    }
}
pub fn event(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
pub fn fixture(derived: bool) -> (Workspace, MaterialOrigin, EventBinding, Molecule, Molecule) {
    let stock = Molecule::dna("Original literal32", STOCK, Topology::Linear).unwrap();
    let circle = Molecule::dna("Independent sealed circle", STOCK, Topology::Circular).unwrap();
    let target = Molecule::dna("Independent edited circle", INTENDED, Topology::Circular).unwrap();
    let w = Workspace::from_documents(
        "Saved inversePCR sources",
        vec![stock.clone(), circle.clone(), target.clone()],
        1,
    )
    .unwrap();
    let p = DesignPlan {
        name: "Actual PCR, phosphorylation and seal".into(),
        target_revision: circle.revision(),
        inputs: vec![PlanInput {
            id: "stock".into(),
            revision: stock.revision(),
            molecule: stock,
        }],
        steps: vec![
            Operation::Pcr {
                template: Source::Input("stock".into()),
                forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
                reverse: Primer::new("R", "TACAGGTC", "").unwrap(),
                template_strand: collagene_core::chemistry::Strand::Top,
            },
            Operation::TreatEnds {
                template: Source::Step(0),
                treatment: EndTreatment::Phosphorylate,
            },
            Operation::Ligate {
                fragments: vec![OrientedSource {
                    source: Source::Step(1),
                    reverse: false,
                }],
                circular: true,
                allow_nicked: false,
                options: Default::default(),
            },
        ],
    };
    let saved = cloning_outputs::finalize(
        &w,
        &p,
        event(&w, 1),
        vec![InputEventBinding {
            input_id: "stock".into(),
            event: event(&w, 0),
        }],
        &cloning_outputs::Options {
            name: "Saved circular template".into(),
            show_prepared: false,
            save_primers: false,
            allow_incomplete: false,
        },
    )
    .unwrap();
    let e = saved.receipt.construct;
    let mut w = saved.workspace;
    let origin = MaterialOrigin::Run(
        run_material::capsule(
            &w,
            &run_material::Selection {
                run_id: saved.receipt.run.run_id,
                step_index: 2,
                event: e.clone(),
            },
        )
        .unwrap(),
    );
    let (e, origin) = if derived {
        w.active_document = e.document_id.clone();
        let catalog = catalog_digest::builtin_catalog().unwrap();
        let id = catalog
            .edition()
            .enzymes
            .iter()
            .find(|e| e.name == "EcoRI")
            .unwrap()
            .id
            .clone();
        let r = derive::Request {
            sources: vec![collagene_core::material_derivation::BoundSource {
                role: collagene_core::material_derivation::SourceRole::Template,
                event: e,
                origin: derive::SourceOrigin::Saved(origin),
            }],
            derivation: derive::Derivation::CatalogDigest {
                selection: catalog_digest::select(&catalog, &[id]).unwrap(),
                options: collagene_core::modified_duplex::checked_inspection_options(),
                fragment_index: 0,
            },
            name: "Saved uncut circle".into(),
            save_primers: false,
        };
        let next = derive::prepare_install(&w, &r).unwrap();
        w = next.workspace;
        let e = next.receipt.product;
        let origin = MaterialOrigin::Derived(derive::capsule(&w, &e).unwrap());
        (e, origin)
    } else {
        (e, origin)
    };
    let source = w.revisions[&e.revision].clone();
    assert_eq!(source.sequence.as_str(), STOCK);
    assert_eq!(source.topology, Topology::Circular);
    (w, origin, e, source, target)
}
