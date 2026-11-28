use collagene_core::{
    documents::workspace::{RevisionEvent, Workspace},
    material::lineage::{self as run_lineage, EventBinding, InputEventBinding},
    output::quantities as q,
    planning::plan::{DesignPlan, Operation, PlanInput, Source},
    sequence::analysis::Primer,
    verification,
};
use collagene_model::{Molecule, Topology};
pub fn pin(w: &Workspace, index: usize) -> EventBinding {
    let d = &w.documents[index];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
pub fn volume(microliters: u64) -> q::units::Volume {
    q::units::Volume::from_unit(
        q::units::Rational::integer(microliters),
        q::units::VolumeUnit::Microliter,
    )
    .unwrap()
}
pub fn fixture() -> (Workspace, q::Request) {
    let source = Molecule::dna("Literal stock", "AGTCCGATGCTACGTA", Topology::Linear).unwrap();
    let target = Molecule::dna(
        "Independent PCR target",
        "AGTCCGATGCTACGTA",
        Topology::Linear,
    )
    .unwrap();
    let private = Molecule::dna(
        "UNRELATED PRIVATE CIRCLE",
        "AACCGGTTAACCGGTT",
        Topology::Circular,
    )
    .unwrap();
    let w = Workspace::from_documents(
        "Explicit arithmetic, not stock authority",
        vec![
            source.clone(),
            target.clone(),
            source.clone(),
            target.clone(),
            private,
        ],
        1,
    )
    .unwrap();
    let plan = DesignPlan {
        name: "Saved native PCR, unchanged status".into(),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "template".into(),
            revision: source.revision(),
            molecule: source,
        }],
        steps: vec![Operation::Pcr {
            template: Source::Input("template".into()),
            forward: Primer::new("F", "AGTCCGAT", "").unwrap(),
            reverse: Primer::new("R", "TACGTAGC", "").unwrap(),
            template_strand: collagene_core::chemistry::Strand::Top,
        }],
    };
    let inputs = vec![InputEventBinding {
        input_id: "template".into(),
        event: pin(&w, 0),
    }];
    let target = pin(&w, 1);
    let saved = run_lineage::record_run(&w, &plan, target.clone(), inputs.clone(), false).unwrap();
    let request = q::Request {
        name: "PRIVATE EXPLICIT QUANTITY PLAN".into(),
        plan_source: q::PlanSource::SavedRun {
            run_id: saved.receipt.run_id,
        },
        plan_revision: collagene_core::planning::plan::DesignPlan::revision(&plan),
        target,
        inputs,
        reactions: vec![q::Reaction {
            id: "pcr".into(),
            step_index: 0,
            final_volume: Some(volume(20)),
            fill: q::FillPolicy::ToFinal {
                component_id: "water".into(),
                label: "Water explicitly declared by operator".into(),
            },
            components: vec![
                q::Component {
                    id: "template-stock".into(),
                    label: "Template user-declared volume".into(),
                    source: Some(Source::Input("template".into())),
                    volume: q::VolumeSpec::Direct {
                        volume: Some(volume(2)),
                    },
                },
                q::Component {
                    id: "other".into(),
                    label: "Other explicitly specified component".into(),
                    source: None,
                    volume: q::VolumeSpec::Direct {
                        volume: Some(volume(3)),
                    },
                },
            ],
        }],
    };
    (saved.workspace, request)
}
#[allow(dead_code)] // Only the native transport consumer advances these shared fixtures.
pub fn advance(w: &mut Workspace, index: usize, name: &str) {
    let m = Molecule::dna(name, "TTTTCCCCAAAAGGGG", Topology::Linear).unwrap();
    let doc = &mut w.documents[index];
    let e = RevisionEvent::new(
        doc.id.clone(),
        m.revision(),
        vec![doc.head_event.clone().unwrap()],
        "Later unrelated head edit",
    )
    .unwrap();
    doc.head = m.revision();
    doc.head_event = Some(e.id.clone());
    w.revisions.insert(m.revision(), m);
    w.history.push(e);
    w.validate().unwrap();
}
