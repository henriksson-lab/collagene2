//! Independent literal activated-vector declaration shared by integration fixtures.
use collagene_core::cloning::topo::physical as p;
use collagene_core::{
    cloning::topo::{self as t, graph_history as g},
    documents::workspace::Workspace,
    material::lineage::EventBinding,
};
use collagene_model::{Molecule, Topology};
pub const FINAL: &str = "AGGGGCATCCCTTCGAATGA";
pub fn event(w: &Workspace, i: usize) -> EventBinding {
    let d = &w.documents[i];
    EventBinding {
        document_id: d.id.clone(),
        event_id: d.head_event.clone().unwrap(),
        revision: d.head.clone(),
    }
}
fn strand(role: &str, length: usize, active: bool) -> t::PreparedStrand {
    t::PreparedStrand {
        id: role.into(),
        segments: vec![t::SourceInterval {
            source_role: role.into(),
            start: 0,
            length,
            reverse: false,
        }],
        five_prime: if active {
            p::FivePrime::Phosphate
        } else {
            p::FivePrime::Hydroxyl
        },
        three_prime: if active {
            p::ThreePrime::Phosphotyrosyl {
                enzyme: p::Enzyme::VacciniaTypeIBTyr274,
                activation_id: format!("activation-{role}"),
            }
        } else {
            p::ThreePrime::Hydroxyl
        },
    }
}
pub fn fixture() -> (Workspace, g::CaptureRequest) {
    let w = Workspace::from_documents(
        "Literal whole graph",
        vec![
            Molecule::dna("vt", "AGGGGCATCCCTT", Topology::Linear).unwrap(),
            Molecule::dna("vb", "AGGGATGCCCCTT", Topology::Linear).unwrap(),
            Molecule::dna("it", "CGAATGA", Topology::Linear).unwrap(),
            Molecule::dna("ib", "CATTCGA", Topology::Linear).unwrap(),
            Molecule::dna("Independent final target", FINAL, Topology::Circular).unwrap(),
            Molecule::dna("UNRELATED PRIVATE", "CCCCAAAAGGGGTTTT", Topology::Circular).unwrap(),
        ],
        0,
    )
    .unwrap();
    let r=g::CaptureRequest{name:"Declared activated vector".into(),bindings:vec![g::SourceBinding{role:"vt".into(),event:event(&w,0),origin:None},g::SourceBinding{role:"vb".into(),event:event(&w,1),origin:None}],request:t::Request{recipe:t::GenerationRecipe{schema_version:1,materials:vec![t::MaterialRecipe::Prepared{declaration_id:"literal activated vector".into(),assumption:"Prepared phosphotyrosyl termini are explicitly declared hypotheses, not evidence of manufacture".into(),strands:vec![strand("vt",13,true),strand("vb",13,true)],pairing:vec![p::PairSpan{strand_a:0,start_a:0,strand_b:1,start_b:0,length:12}]}]},steps:vec![],engine_limits:p::Limits::default(),projection_options:p::projection::Options::default()}};
    (w, r)
}
pub fn captured() -> g::PreparedCapture {
    let (w, r) = fixture();
    g::prepare_capture(&w, &r).unwrap()
}
pub fn partial(w: &Workspace, selection: &g::Selection) -> g::CaptureRequest {
    g::CaptureRequest {
        name: "One transfer, entire remaining graph".into(),
        bindings: vec![
            g::SourceBinding {
                role: "it".into(),
                event: event(w, 2),
                origin: None,
            },
            g::SourceBinding {
                role: "ib".into(),
                event: event(w, 3),
                origin: None,
            },
        ],
        request: t::Request {
            recipe: t::GenerationRecipe {
                schema_version: 1,
                materials: vec![
                    t::MaterialRecipe::SavedGraph {
                        material_id: "retained vector".into(),
                        origin: g::documentary_capsule(w, selection).unwrap(),
                    },
                    t::MaterialRecipe::Prepared {
                        declaration_id: "literal insert".into(),
                        assumption:
                            "Explicit OH insert strands, not activation inferred from sequence"
                                .into(),
                        strands: vec![strand("it", 7, false), strand("ib", 7, false)],
                        pairing: vec![p::PairSpan {
                            strand_a: 0,
                            start_a: 0,
                            strand_b: 1,
                            start_b: 0,
                            length: 6,
                        }],
                    },
                ],
            },
            steps: vec![p::Step::Transfer(p::Transfer {
                donor: p::Call {
                    strand: 0,
                    index: 12,
                },
                acceptor: p::Call {
                    strand: 2,
                    index: 0,
                },
                geometry: p::Geometry::ThreePrimeTA,
            })],
            engine_limits: p::Limits::default(),
            projection_options: p::projection::Options::default(),
        },
    }
}
pub fn continuation(w: &Workspace, selection: &g::Selection) -> t::Request {
    t::Request {
        recipe: t::GenerationRecipe {
            schema_version: 1,
            materials: vec![t::MaterialRecipe::SavedGraph {
                material_id: "entire partial graph".into(),
                origin: g::documentary_capsule(w, selection).unwrap(),
            }],
        },
        steps: vec![p::Step::Transfer(p::Transfer {
            donor: p::Call {
                strand: 1,
                index: 12,
            },
            acceptor: p::Call {
                strand: 3,
                index: 0,
            },
            geometry: p::Geometry::ThreePrimeTA,
        })],
        engine_limits: p::Limits::default(),
        projection_options: p::projection::Options::default(),
    }
}
