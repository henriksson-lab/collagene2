//! Independent old positional DTO encoders: fields after each nested scientific
//! report carry sentinels, so an EOF/default shortcut cannot satisfy these tests.
pub use collagene_core::*;
use serde::{Deserialize, Serialize};
use serde_json::json;
#[path = "../src/protocol_legacy_v64.rs"]
mod decoder;
#[derive(Clone, Serialize, Deserialize)]
struct OldReport {
    plan_revision: String,
    target_revision: String,
    engine_version: String,
    status: verification::Status,
    findings: Vec<String>,
    assumptions: Vec<String>,
    products: Vec<Molecule>,
    exhaustive_alternative_search: bool,
    step_reports: Vec<verification::StepReport>,
    annotation_reports: Vec<verification::AnnotationStepReport>,
}
macro_rules! old {($name:ident<R>{$($field:ident:$ty:ty),*$(,)?})=>{#[derive(Serialize,Deserialize)]struct $name<R=OldReport>{$($field:$ty),*}};($name:ident{$($field:ident:$ty:ty),*$(,)?})=>{#[derive(Serialize,Deserialize)]struct $name{$($field:$ty),*}};}
#[derive(Serialize, Deserialize)]
struct OldInverse<R = OldReport> {
    source_revision: String,
    target_revision: String,
    options: mutagenesis::InversePcrOptions,
    forward: analysis::Primer,
    reverse: analysis::Primer,
    forward_location: Location,
    reverse_location: Location,
    retained_template_length: usize,
    plan: verification::DesignPlan,
    verification: R,
    warnings: Vec<String>,
}
#[derive(Serialize, Deserialize)]
struct OldReplay<R = OldReport> {
    run_id: String,
    source_workspace_fingerprint: String,
    recorded_status: verification::Status,
    recorded_engine: String,
    current_report: R,
    engine_matches: bool,
    products_match: bool,
    status_matches: bool,
}
old!(OldOverlap<R>{context_revision:String,compiled:assembly_recipe::Compiled,verification:R,prepared:Vec<Molecule>,profile:Option<overlap_profiles::ProfileReport>,preparation_complete:bool});
old!(OldLigation<R>{context_revision:String,compiled:ligation_recipe::Compiled,verification:R,prepared:Vec<Molecule>,physical:Vec<Option<chemistry::DuplexSummary>>,product_physical:Option<chemistry::DuplexSummary>,junctions:Vec<ligation_recipe::Junction>,preparation_complete:bool});
old!(OldJunction<R>{engine_version:String,context_revision:String,plan_revision:String,target_revision:String,options:assembly_junctions::Options,verification:R,product_revision:Option<String>,product_length:Option<usize>,product_topology:Option<Topology>,method:Option<assembly_junctions::Method>,junctions:Vec<assembly_junctions::Junction>,complete:bool,omitted_junctions:usize,warnings:Vec<String>});
old!(OldPool<R>{context_revision:String,request:prepared_golden_gate::Request,target:Molecule,preview:OldLigation<R>,missing_physical:Vec<usize>,pool:Option<golden_gate::GoldenGateReport>,assessments:Vec<prepared_golden_gate::Assessment>,warnings:Vec<String>});
old!(OldReaction<R>{context_revision:String,request:golden_gate_reaction::Request,target:Molecule,probe_plan:verification::DesignPlan,materials:Vec<golden_gate_reaction::Material>,verification:R,pool:Option<golden_gate::GoldenGateReport>,assessments:Vec<prepared_golden_gate::Assessment>,warnings:Vec<String>});
old!(OldDomestication<R>{context_revision:String,request:domestication_preparation::Request,source_revision:String,target:Molecule,check:domestication::CheckReport,preparation:mutagenesis::InversePcrPreparation,verification:R,warnings:Vec<String>});
#[derive(Serialize)]
struct Wire<'a, T> {
    id: u32,
    revision: &'a str,
    tag: u32,
    payload: T,
}
fn envelope(payload: &[u8], version: u16) -> Vec<u8> {
    let mut out = b"CGR\0".to_vec();
    out.extend_from_slice(&version.to_le_bytes());
    out.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    out.extend_from_slice(blake3::hash(payload).as_bytes());
    out.extend_from_slice(payload);
    out
}
fn encode<T: Serialize>(tag: u32, payload: &T, version: u16) -> Vec<u8> {
    envelope(
        &postcard::to_allocvec(&Wire {
            id: 733,
            revision: "retained-old-context",
            tag,
            payload,
        })
        .unwrap(),
        version,
    )
}
fn report() -> OldReport {
    OldReport {
        plan_revision: "plan-source-sentinel".into(),
        target_revision: "independent-target-sentinel".into(),
        engine_version: "0.1.0/science-15".into(),
        status: verification::Status::Incomplete,
        findings: vec!["old findings retained".into()],
        assumptions: vec!["old assumptions retained".into()],
        products: vec![molecule()],
        exhaustive_alternative_search: false,
        step_reports: vec![verification::StepReport {
            step_index: 3,
            topo: None,
        }],
        annotation_reports: vec![],
    }
}
fn molecule() -> Molecule {
    Molecule::dna("retained original document", "ACGTAG", Topology::Circular).unwrap()
}
fn plan() -> verification::DesignPlan {
    verification::DesignPlan {
        name: "old positional plan".into(),
        target_revision: molecule().revision(),
        inputs: vec![],
        steps: vec![],
    }
}
fn ligation() -> OldLigation {
    OldLigation {
        context_revision: "ligation-context".into(),
        compiled: ligation_recipe::Compiled {
            plan: plan(),
            rows: vec![],
            ligation_step: 9,
        },
        verification: report(),
        prepared: vec![molecule()],
        physical: vec![None],
        product_physical: None,
        junctions: vec![ligation_recipe::Junction {
            index: 7,
            upstream: 2,
            downstream: 3,
            closing: true,
            upstream_end: None,
            downstream_end: None,
        }],
        preparation_complete: true,
    }
}
fn compare<T: Serialize>(tag: u32, value: T, path: &str) {
    let original = serde_json::to_value(&value).unwrap();
    let bytes = encode(tag, &value, 64);
    let response = decoder::decode(&bytes, 32 << 20).unwrap().unwrap();
    assert_eq!(response.id, 733);
    assert_eq!(response.revision, "retained-old-context");
    let relabelled = envelope(&postcard::to_allocvec(&response).unwrap(), 64);
    assert!(
        decoder::decode(&relabelled, 32 << 20).is_err(),
        "new report field cannot be decoded as old tag {tag}"
    );
    let actual = serde_json::to_value(&response.result).unwrap();
    let actual = actual.as_object().unwrap().values().next().unwrap();
    let mut expected = original.clone();
    let report = if path.is_empty() {
        &mut expected
    } else {
        expected.pointer_mut(path).unwrap()
    };
    report
        .as_object_mut()
        .unwrap()
        .insert("processed_step_reports".into(), json!([]));
    report
        .as_object_mut()
        .unwrap()
        .insert("physical_processing_reports".into(), json!([]));
    assert_eq!(*actual, expected, "legacy result tag {tag}");
    let mut trailing = postcard::to_allocvec(&Wire {
        id: 733,
        revision: "retained-old-context",
        tag,
        payload: value,
    })
    .unwrap();
    trailing.push(0);
    assert!(
        decoder::decode(&envelope(&trailing, 64), 32 << 20).is_err(),
        "trailing fields tag {tag}"
    );
    upgrade65(tag, original, path);
}
#[test]
fn bare_legacy_report_preserves_engine_source_product_and_empty_new_evidence() {
    compare(4, report(), "");
    let migrated = protocol::decode_response(&encode(4, &report(), 64)).unwrap();
    let protocol::AnalysisResult::Verification(migrated) = migrated.result else {
        panic!("expected legacy verification dispatch")
    };
    assert_eq!(migrated.engine_version, "0.1.0/science-15");
    assert!(migrated.processed_step_reports.is_empty());
}
#[test]
fn old_inverse_pcr_preserves_following_warning_fields() {
    let primer = analysis::Primer {
        name: "old primer".into(),
        annealing: collagene_core::Sequence::dna("ACGT").unwrap(),
        tail: "GG".into(),
    };
    let location = Location::range(0, 4);
    compare(
        16,
        OldInverse {
            source_revision: "source".into(),
            target_revision: "target".into(),
            options: mutagenesis::InversePcrOptions::default(),
            forward: primer.clone(),
            reverse: primer,
            forward_location: location.clone(),
            reverse_location: location,
            retained_template_length: 17,
            plan: plan(),
            verification: report(),
            warnings: vec!["after old inverse report".into()],
        },
        "/verification",
    );
}
#[test]
fn old_replay_preserves_three_following_booleans() {
    compare(
        41,
        OldReplay {
            run_id: "old run".into(),
            source_workspace_fingerprint: "original workspace".into(),
            recorded_status: verification::Status::Incomplete,
            recorded_engine: "science-14".into(),
            current_report: report(),
            engine_matches: false,
            products_match: true,
            status_matches: false,
        },
        "/current_report",
    );
}
#[test]
fn old_overlap_preview_preserves_documents_and_completion_after_report() {
    compare(
        86,
        OldOverlap {
            context_revision: "overlap-context".into(),
            compiled: assembly_recipe::Compiled {
                plan: plan(),
                rows: vec![],
                assembly_step: 5,
            },
            verification: report(),
            prepared: vec![molecule()],
            profile: None,
            preparation_complete: true,
        },
        "/verification",
    );
}
#[test]
fn old_ligation_preview_preserves_physical_vector_and_junctions() {
    compare(87, ligation(), "/verification");
}
#[test]
fn old_junction_inspection_preserves_all_post_report_fields() {
    compare(
        88,
        OldJunction {
            engine_version: "assembly-junctions-1".into(),
            context_revision: "junction-context".into(),
            plan_revision: "plan".into(),
            target_revision: "target".into(),
            options: assembly_junctions::Options::default(),
            verification: report(),
            product_revision: Some("actual-product".into()),
            product_length: Some(987),
            product_topology: Some(Topology::Circular),
            method: Some(assembly_junctions::Method::Ligate),
            junctions: vec![],
            complete: true,
            omitted_junctions: 3,
            warnings: vec!["after junction report".into()],
        },
        "/verification",
    );
}
#[test]
fn old_prepared_pool_migrates_double_nested_ligation_report() {
    compare(
        105,
        OldPool {
            context_revision: "pool-context".into(),
            request: prepared_golden_gate::Request {
                recipe: ligation_recipe::Recipe {
                    name: "pool recipe".into(),
                    fragments: vec![],
                    circular: true,
                    allow_nicked: false,
                },
                options: golden_gate::GoldenGateOptions::default(),
            },
            target: molecule(),
            preview: ligation(),
            missing_physical: vec![7, 9],
            pool: None,
            assessments: vec![],
            warnings: vec!["after nested pool report".into()],
        },
        "/preview/verification",
    );
}
#[test]
fn old_reaction_preserves_following_assessments_and_warning() {
    compare(
        107,
        OldReaction {
            context_revision: "reaction-context".into(),
            request: golden_gate_reaction::Request {
                name: "reaction".into(),
                rows: vec![],
                options: golden_gate::GoldenGateOptions::default(),
            },
            target: molecule(),
            probe_plan: plan(),
            materials: vec![],
            verification: report(),
            pool: None,
            assessments: vec![prepared_golden_gate::Assessment {
                candidate_index: 3,
                directed_sequence_key: "directed-key".into(),
                physical_sequence_key: "physical-key".into(),
                target_matches: false,
                source_roles: vec![2, 4],
            }],
            warnings: vec!["after reaction report".into()],
        },
        "/verification",
    );
}
#[test]
fn old_domestication_preparation_preserves_independent_check_and_final_warnings() {
    let source = molecule();
    let design = domestication::Request {
        source_revision: source.revision(),
        selection: catalog_digest::select(
            &catalog_digest::builtin_catalog().unwrap(),
            &["bsai".into()],
        )
        .unwrap(),
        scope: vec![domestication::SourceRange { start: 0, end: 6 }],
        retained_sites: vec![],
        coding_features: vec![],
        protected: vec![],
        options: codon_design::DesignOptions::default(),
        limits: domestication::Limits::default(),
    };
    let selected = domestication::Selection {
        context_revision: "old proposed context".into(),
        candidate_id: "old literal candidate".into(),
        patches: vec![],
    };
    let check = domestication::CheckReport {
        context_revision: "old checked context".into(),
        request: design.clone(),
        selection: selected.clone(),
        source_revision: source.revision(),
        target_revision: "original independent intent".into(),
        candidate_id: selected.candidate_id.clone(),
        sites_before: vec![],
        sites_after: vec![],
        features: vec![],
        warnings: vec!["original check annotation".into()],
    };
    let primer = analysis::Primer {
        name: "literal orderable primer".into(),
        annealing: collagene_core::Sequence::dna("ACGT").unwrap(),
        tail: "C".into(),
    };
    let preparation = mutagenesis::InversePcrPreparation {
        source_revision: source.revision(),
        target_revision: "original independent intent".into(),
        options: mutagenesis::InversePcrOptions::default(),
        forward: primer.clone(),
        reverse: primer,
        forward_location: Location::range(0, 4),
        reverse_location: Location::range(2, 6),
        retained_template_length: 4,
        plan: plan(),
        warnings: vec!["literal preparation warning".into()],
    };
    compare(
        134,
        OldDomestication {
            context_revision: "physical-context".into(),
            request: domestication_preparation::Request {
                design,
                selection: selected,
                route: domestication_preparation::Route::InversePcr {
                    edit_start: 2,
                    edit_end: 3,
                    forward_length: 4,
                    reverse_length: 4,
                    reverse_tail_bases: 0,
                },
            },
            source_revision: source.revision(),
            target: source,
            check,
            preparation,
            verification: report(),
            warnings: vec!["AFTER old domestication report".into()],
        },
        "/verification",
    );
}
#[test]
fn envelopes_are_checked_before_payloads_and_unrelated_results_fall_through() {
    let valid = encode(4, &report(), 64);
    assert!(decoder::decode(&valid, 1).is_err());
    let mut wrong = valid.clone();
    wrong[0] = b'X';
    assert!(decoder::decode(&wrong, 32 << 20).is_err());
    let mut wrong = valid.clone();
    wrong[10] ^= 1;
    assert!(decoder::decode(&wrong, 32 << 20).is_err());
    let mut wrong = valid.clone();
    wrong.push(0);
    assert!(decoder::decode(&wrong, 32 << 20).is_err());
    let mut wrong = valid.clone();
    wrong[6..10].copy_from_slice(&u32::MAX.to_le_bytes());
    assert!(decoder::decode(&wrong, 32 << 20).is_err());
    let mut wrong = valid.clone();
    wrong[4..6].copy_from_slice(&0u16.to_le_bytes());
    assert!(decoder::decode(&wrong, 32 << 20).is_err());
    let mut wrong = valid.clone();
    wrong[4..6].copy_from_slice(&(protocol::WORKER_SCHEMA_VERSION + 1).to_le_bytes());
    assert!(decoder::decode(&wrong, 32 << 20).is_err());
    // A known unrelated historical Error is left for the ordinary protocol decoder.
    assert!(
        decoder::decode(&encode(2, &"historical error", 64), 32 << 20)
            .unwrap()
            .is_none()
    );
    // Earliest bare layouts belong to their own existing legacy migrations.
    assert!(
        decoder::decode(&encode(4, &report(), 26), 32 << 20)
            .unwrap()
            .is_none()
    );
    // Current envelopes are decoded by the current schema implementation.
    assert!(
        decoder::decode(
            &encode(4, &report(), protocol::WORKER_SCHEMA_VERSION),
            32 << 20
        )
        .unwrap()
        .is_none()
    );
}
// Actual nested report layouts used before graph evidence (eight fields) and
// before annotations (nine fields). These encoders do not serialize a current
// report under an old header.
old!(EightFieldReport{plan_revision:String,target_revision:String,engine_version:String,status:verification::Status,findings:Vec<String>,assumptions:Vec<String>,products:Vec<Molecule>,exhaustive_alternative_search:bool});
old!(NineFieldReport{plan_revision:String,target_revision:String,engine_version:String,status:verification::Status,findings:Vec<String>,assumptions:Vec<String>,products:Vec<Molecule>,exhaustive_alternative_search:bool,step_reports:Vec<verification::StepReport>});
fn old_era_report<R: for<'de> Deserialize<'de>>(version: u16) -> R {
    let mut value = serde_json::to_value(report()).unwrap();
    let object = value.as_object_mut().unwrap();
    object.insert("plan_revision".into(), json!("1".repeat(64)));
    object.insert("target_revision".into(), json!("2".repeat(64)));
    object.insert(
        "engine_version".into(),
        json!(if version < 25 {
            "0.1.0/science-9"
        } else {
            "0.1.0/science-10"
        }),
    );
    object.remove("annotation_reports");
    if version < 25 {
        object.remove("step_reports");
    }
    serde_json::from_value(value).unwrap()
}
fn public_era_roundtrip<T: Serialize>(version: u16, tag: u32, payload: T, path: &str) {
    let mut expected = serde_json::to_value(&payload).unwrap();
    let p = expected.pointer_mut(path).unwrap().as_object_mut().unwrap();
    p.entry("step_reports").or_insert(json!([]));
    p.insert("annotation_reports".into(), json!([]));
    p.insert("processed_step_reports".into(), json!([]));
    p.insert("physical_processing_reports".into(), json!([]));
    let bytes = envelope(
        &postcard::to_allocvec(&Wire {
            id: 912,
            revision: &"b".repeat(64),
            tag,
            payload,
        })
        .unwrap(),
        version,
    );
    let actual = protocol::decode_response(&bytes).unwrap();
    assert_eq!(actual.id, 912);
    assert_eq!(actual.revision, "b".repeat(64));
    let value = serde_json::to_value(actual.result).unwrap();
    assert_eq!(
        value.as_object().unwrap().values().next().unwrap(),
        &expected,
        "old nested layout tag{tag} schema{version}"
    );
}
fn inverse_era<R>(report: R) -> OldInverse<R> {
    let primer = analysis::Primer {
        name: "historical primer".into(),
        annealing: collagene_core::Sequence::dna("ACGT").unwrap(),
        tail: "GG".into(),
    };
    OldInverse {
        source_revision: "source-era".into(),
        target_revision: "independent-target-era".into(),
        options: mutagenesis::InversePcrOptions::default(),
        forward: primer.clone(),
        reverse: primer,
        forward_location: Location::range(0, 4),
        reverse_location: Location::range(2, 6),
        retained_template_length: 12,
        plan: plan(),
        verification: report,
        warnings: vec!["TAIL WARNING AFTER EIGHT/NINE-FIELD REPORT".into()],
    }
}
fn replay_era<R>(report: R) -> OldReplay<R> {
    OldReplay {
        run_id: "a".repeat(64),
        source_workspace_fingerprint: "b".repeat(64),
        recorded_status: verification::Status::ConsistentUnderAssumptions,
        recorded_engine: "0.1.0/science-8".into(),
        current_report: report,
        engine_matches: false,
        products_match: true,
        status_matches: false,
    }
}
#[test]
fn historical_inverse_nested_eight_and_nine_field_reports_publicly_migrate() {
    public_era_roundtrip(
        24,
        16,
        inverse_era(old_era_report::<EightFieldReport>(24)),
        "/verification",
    );
    for version in [25, 26] {
        public_era_roundtrip(
            version,
            16,
            inverse_era(old_era_report::<NineFieldReport>(version)),
            "/verification",
        );
    }
}
#[test]
fn historical_replay_nested_eight_and_nine_field_reports_preserve_recorded_status_and_tail_flags() {
    public_era_roundtrip(
        24,
        41,
        replay_era(old_era_report::<EightFieldReport>(24)),
        "/current_report",
    );
    for version in [25, 26] {
        public_era_roundtrip(
            version,
            41,
            replay_era(old_era_report::<NineFieldReport>(version)),
            "/current_report",
        );
    }
}
old!(ElevenFieldReport{plan_revision:String,target_revision:String,engine_version:String,status:verification::Status,findings:Vec<String>,assumptions:Vec<String>,products:Vec<Molecule>,exhaustive_alternative_search:bool,step_reports:Vec<verification::StepReport>,annotation_reports:Vec<verification::AnnotationStepReport>,processed_step_reports:Vec<verification::ProcessedTopoStepReport>});
fn real_old_processed_evidence() -> topo_processed_product::Report {
    use topo::physical as p;
    let dna = ["AGGGGCATCCCTT", "AGGGATGCCCCTT", "CGAATGA", "CATTCGA"];
    let inputs = dna
        .iter()
        .enumerate()
        .map(|(i, dna)| topo::Input {
            role: format!("old-source{i}"),
            molecule: Molecule::dna(format!("original{i}"), dna, Topology::Linear).unwrap(),
        })
        .collect::<Vec<_>>();
    let strands = dna
        .iter()
        .enumerate()
        .map(|(i, dna)| topo::PreparedStrand {
            id: format!("strand{i}"),
            segments: vec![topo::SourceInterval {
                source_role: format!("old-source{i}"),
                start: 0,
                length: dna.len(),
                reverse: false,
            }],
            five_prime: if i < 2 {
                p::FivePrime::Phosphate
            } else {
                p::FivePrime::Hydroxyl
            },
            three_prime: if i < 2 {
                p::ThreePrime::Phosphotyrosyl {
                    enzyme: p::Enzyme::VacciniaTypeIBTyr274,
                    activation_id: format!("activation{i}"),
                }
            } else {
                p::ThreePrime::Hydroxyl
            },
        })
        .collect();
    let preparation = topo::Request {
        recipe: topo::GenerationRecipe {
            schema_version: 1,
            materials: vec![topo::MaterialRecipe::Prepared {
                declaration_id: "historical explicit input geometry".into(),
                assumption: "Declared original activation hypothesis".into(),
                strands,
                pairing: vec![
                    p::PairSpan {
                        strand_a: 0,
                        start_a: 0,
                        strand_b: 1,
                        start_b: 0,
                        length: 12,
                    },
                    p::PairSpan {
                        strand_a: 2,
                        start_a: 0,
                        strand_b: 3,
                        start_b: 0,
                        length: 6,
                    },
                ],
            }],
        },
        transfers: vec![
            p::Transfer {
                donor: 0,
                acceptor: 2,
                geometry: p::Geometry::ThreePrimeTA,
            },
            p::Transfer {
                donor: 1,
                acceptor: 3,
                geometry: p::Geometry::ThreePrimeTA,
            },
        ],
        engine_limits: p::Limits::default(),
        projection_options: p::projection::Options::default(),
    };
    let sealing = p::nick_sealing::Request {
        ligase_assumption: "Historical ideal ligase hypothesis".into(),
        nicks: vec![
            p::nick_sealing::NickSelection {
                three: p::Call {
                    strand: 2,
                    index: 6,
                },
                five: p::Call {
                    strand: 0,
                    index: 0,
                },
            },
            p::nick_sealing::NickSelection {
                three: p::Call {
                    strand: 3,
                    index: 6,
                },
                five: p::Call {
                    strand: 1,
                    index: 0,
                },
            },
        ],
    };
    let prepared = topo::simulate(&inputs, &preparation).unwrap().preparation;
    let context = p::nick_sealing::revision(
        &prepared.substrates,
        &preparation.transfers,
        &preparation.engine_limits,
        &sealing,
    )
    .unwrap();
    topo_processed_product::inspect(
        &inputs,
        &topo_processed_product::Request {
            treatment: topo_nick_sealing::Request {
                preparation,
                sealing,
            },
            selection: p::processed_projection::Selection {
                sealing_context_revision: context,
                component: 0,
                reference_seed: p::Call {
                    strand: 0,
                    index: 0,
                },
            },
            options: p::projection::Options::default(),
        },
    )
    .unwrap()
}
fn eleven_field_report() -> ElevenFieldReport {
    let mut value = serde_json::to_value(report()).unwrap();
    value["engine_version"] = json!("0.1.0/science-16");
    value["processed_step_reports"] = json!([verification::ProcessedTopoStepReport {
        step_index: 0,
        evidence: real_old_processed_evidence()
    }]);
    serde_json::from_value(value).unwrap()
}
fn compare65<T: Serialize>(tag: u32, value: T, path: &str) {
    let mut expected = serde_json::to_value(&value).unwrap();
    let old = if path.is_empty() {
        &mut expected
    } else {
        expected.pointer_mut(path).unwrap()
    };
    assert_eq!(old["processed_step_reports"].as_array().unwrap().len(), 1);
    old.as_object_mut()
        .unwrap()
        .insert("physical_processing_reports".into(), json!([]));
    let oldbytes = encode(tag, &value, 65);
    let response = decoder::decode(&oldbytes, 32 << 20).unwrap().unwrap();
    let actual = serde_json::to_value(&response.result).unwrap();
    assert_eq!(
        actual.as_object().unwrap().values().next().unwrap(),
        &expected
    );
    let current_as_old = envelope(&postcard::to_allocvec(&response).unwrap(), 65);
    assert!(decoder::decode(&current_as_old, 32 << 20).is_err());
    let mut trailing = postcard::to_allocvec(&Wire {
        id: 733,
        revision: "retained-old-context",
        tag,
        payload: value,
    })
    .unwrap();
    trailing.push(0);
    assert!(decoder::decode(&envelope(&trailing, 65), 32 << 20).is_err());
}
#[test]
fn actual_schema65_bare_processed_graph_survives_without_invented_generic_actions() {
    let old = eleven_field_report();
    let wire = encode(4, &old, 65);
    let restored = protocol::decode_response(&wire).unwrap();
    let protocol::AnalysisResult::Verification(r) = restored.result else {
        panic!("verification")
    };
    assert_eq!(r.engine_version, "0.1.0/science-16");
    assert_eq!(r.processed_step_reports.len(), 1);
    assert_eq!(
        r.processed_step_reports[0]
            .evidence
            .projection
            .sealing
            .source_bases,
        40
    );
    assert!(r.physical_processing_reports.is_empty());
    compare65(4, old, "");
}
#[test]
fn actual_schema65_nested_processed_graph_preserves_inverse_tail_warnings() {
    compare65(16, inverse_era(eleven_field_report()), "/verification");
}
// Each independently encoded <=64 outer fixture also traverses its real65
// report era with nonempty original processed graph evidence and trailing fields.
fn upgrade65(tag: u32, mut value: serde_json::Value, path: &str) {
    let report = serde_json::to_value(eleven_field_report()).unwrap();
    if path.is_empty() {
        value = report
    } else {
        *value.pointer_mut(path).unwrap() = report;
    }
    macro_rules! run {
        ($t:ty) => {
            compare65(tag, serde_json::from_value::<$t>(value).unwrap(), path)
        };
    }
    match tag {
        4 => run!(ElevenFieldReport),
        16 => run!(OldInverse<ElevenFieldReport>),
        41 => run!(OldReplay<ElevenFieldReport>),
        86 => run!(OldOverlap<ElevenFieldReport>),
        87 => run!(OldLigation<ElevenFieldReport>),
        88 => run!(OldJunction<ElevenFieldReport>),
        105 => run!(OldPool<ElevenFieldReport>),
        107 => run!(OldReaction<ElevenFieldReport>),
        134 => run!(OldDomestication<ElevenFieldReport>),
        _ => panic!("unknown old layout fixture"),
    }
}
