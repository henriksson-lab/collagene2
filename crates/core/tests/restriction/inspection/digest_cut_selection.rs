use collagene_core::{
    chemistry::{FivePrime, ThreePrime},
    cloning::restriction::{
        catalog as c,
        cut_selection::{self as s, Boundary as B, Choice},
    },
};
use collagene_model::{Molecule, Topology};
fn source(dna: &str, topology: Topology) -> Molecule {
    Molecule::dna("unchanged original", dna, topology).unwrap()
}
fn select(id: &str) -> (c::catalog::Catalog, c::CatalogSelection) {
    let catalog = c::builtin_catalog().unwrap();
    let selection = c::select(&catalog, &[id.into()]).unwrap();
    (catalog, selection)
}
fn cut(site_index: usize) -> B {
    B::Cut { site_index }
}
fn preview(
    m: &Molecule,
    selection: &c::CatalogSelection,
    catalog: &c::catalog::Catalog,
    choice: &Choice,
) -> s::Report {
    let o = collagene_core::cloning::restriction::inspection::InspectionOptions::default();
    let batch = c::batch_revision(catalog, selection, "actual-source-role", m, &o).unwrap();
    s::preview(
        catalog,
        selection,
        "actual-source-role",
        m,
        &o,
        &batch,
        choice,
    )
    .unwrap()
}
#[test]
fn circle_two_actual_cut_instances_select_forward_or_complement_fragment() {
    let (catalog, selection) = select("ecori");
    let m = source("GAATTCAAAAGAATTCCCC", Topology::Circular);
    let before = serde_json::to_value(&m).unwrap();
    let chosen = preview(
        &m,
        &selection,
        &catalog,
        &Choice {
            first: cut(0),
            second: cut(1),
            complement: false,
        },
    );
    assert_eq!(chosen.product.product.fragment_index, 0);
    assert_eq!(
        chosen.product.product.molecule.sequence.as_str(),
        "AATTCAAAAG"
    );
    let f = &chosen.product.product.summary;
    assert_eq!(
        (
            f.top_start,
            f.bottom_start,
            f.duplex.top_length,
            f.duplex.bottom_length
        ),
        (1, 5, 10, 10)
    );
    assert_eq!(f.duplex.bottom_start, 4);
    let ends = f.duplex.ends.as_ref().unwrap();
    assert_eq!(
        (ends.0.overhang.as_str(), ends.1.overhang.as_str()),
        ("AATT", "AATT")
    );
    assert_eq!(ends.0.chemistry.five_prime, FivePrime::Phosphate);
    assert_eq!(ends.1.chemistry.three_prime, ThreePrime::Hydroxyl);
    let other = preview(
        &m,
        &selection,
        &catalog,
        &Choice {
            first: cut(0),
            second: cut(1),
            complement: true,
        },
    );
    assert_eq!(other.product.product.fragment_index, 1);
    assert_eq!(
        other.product.product.molecule.sequence.as_str(),
        "AATTCCCCG"
    );
    assert_eq!(
        (
            other.product.product.summary.top_start,
            other.product.product.summary.bottom_start
        ),
        (11, 15)
    );
    assert_eq!(serde_json::to_value(&m).unwrap(), before);
    assert_ne!(chosen.context_revision, other.context_revision);
}
#[test]
fn one_circle_cut_selects_a_full_turn_with_global_bottom_geometry() {
    let (catalog, selection) = select("nhei");
    let m = source("AAAAGCTAGCCCC", Topology::Circular);
    let chosen = preview(
        &m,
        &selection,
        &catalog,
        &Choice {
            first: cut(0),
            second: cut(0),
            complement: false,
        },
    );
    assert_eq!(
        chosen.product.product.molecule.sequence.as_str(),
        "CTAGCCCCAAAAG"
    );
    let f = &chosen.product.product.summary;
    assert_eq!(
        (
            f.top_start,
            f.bottom_start,
            f.duplex.top_length,
            f.duplex.bottom_length
        ),
        (5, 9, 13, 13)
    );
    assert_eq!(f.duplex.bottom_start, 4);
    let reversedchoice = preview(
        &m,
        &selection,
        &catalog,
        &Choice {
            first: cut(0),
            second: cut(0),
            complement: true,
        },
    );
    assert_eq!(
        chosen.product.product.molecule.sequence,
        reversedchoice.product.product.molecule.sequence
    );
}
#[test]
fn linear_end_fragments_normalize_click_order_and_retain_actual_unequal_strands() {
    let (catalog, selection) = select("ecori");
    let m = source("AAAAGAATTCCCC", Topology::Linear);
    let left = preview(
        &m,
        &selection,
        &catalog,
        &Choice {
            first: B::SourceStart,
            second: cut(0),
            complement: false,
        },
    );
    assert_eq!(left.product.product.molecule.sequence.as_str(), "AAAAG");
    assert_eq!(
        (
            left.product.product.summary.top_start,
            left.product.product.summary.bottom_start,
            left.product.product.summary.duplex.top_length,
            left.product.product.summary.duplex.bottom_length
        ),
        (0, 0, 5, 9)
    );
    let right = preview(
        &m,
        &selection,
        &catalog,
        &Choice {
            first: B::SourceEnd,
            second: cut(0),
            complement: false,
        },
    );
    assert_eq!(right.product.product.molecule.sequence.as_str(), "AATTCCCC");
    assert_eq!(
        (
            right.product.product.summary.top_start,
            right.product.product.summary.bottom_start,
            right.product.product.summary.duplex.top_length,
            right.product.product.summary.duplex.bottom_length
        ),
        (5, 9, 8, 4)
    );
    let options = collagene_core::cloning::restriction::inspection::InspectionOptions::default();
    let batch =
        c::batch_revision(&catalog, &selection, "actual-source-role", &m, &options).unwrap();
    assert!(
        s::preview(
            &catalog,
            &selection,
            "actual-source-role",
            &m,
            &options,
            &batch,
            &Choice {
                first: B::SourceStart,
                second: cut(0),
                complement: true
            }
        )
        .is_err()
    );
}
#[test]
fn internal_extra_cuts_possible_sites_and_stale_batches_cannot_be_sliced_into_a_product() {
    let (catalog, selection) = select("ecori");
    let m = source("GAATTCAAAAGAATTCCCCGAATTCTTT", Topology::Circular);
    let options = collagene_core::cloning::restriction::inspection::InspectionOptions::default();
    let batch =
        c::batch_revision(&catalog, &selection, "actual-source-role", &m, &options).unwrap();
    let choice = Choice {
        first: cut(0),
        second: cut(2),
        complement: false,
    };
    let error = s::preview(
        &catalog,
        &selection,
        "actual-source-role",
        &m,
        &options,
        &batch,
        &choice,
    )
    .unwrap_err();
    assert!(error.to_string().contains("additional internal cuts"));
    assert!(
        s::preview(
            &catalog,
            &selection,
            "actual-source-role",
            &m,
            &options,
            &"0".repeat(64),
            &choice
        )
        .is_err()
    );
    let possible = source("GANTTCAAAA", Topology::Circular);
    let batch = c::batch_revision(
        &catalog,
        &selection,
        "actual-source-role",
        &possible,
        &options,
    )
    .unwrap();
    assert!(
        s::preview(
            &catalog,
            &selection,
            "actual-source-role",
            &possible,
            &options,
            &batch,
            &Choice {
                first: cut(0),
                second: cut(0),
                complement: false
            }
        )
        .is_err()
    );
    let (catalog, selection) = select("bsai");
    let outside = source("GGTCTC", Topology::Linear);
    let batch = c::batch_revision(
        &catalog,
        &selection,
        "actual-source-role",
        &outside,
        &options,
    )
    .unwrap();
    assert!(
        s::preview(
            &catalog,
            &selection,
            "actual-source-role",
            &outside,
            &options,
            &batch,
            &Choice {
                first: B::SourceStart,
                second: cut(0),
                complement: false
            }
        )
        .is_err()
    );
}
#[test]
fn decoded_receipt_choice_and_product_must_remain_bound_without_restoring_state() {
    let (catalog, selection) = select("ecori");
    let m = source("GAATTCAAAAGAATTCCCC", Topology::Circular);
    let choice = Choice {
        first: cut(0),
        second: cut(1),
        complement: false,
    };
    let report = preview(&m, &selection, &catalog, &choice);
    let options = collagene_core::cloning::restriction::inspection::InspectionOptions::default();
    let batch =
        c::batch_revision(&catalog, &selection, "actual-source-role", &m, &options).unwrap();
    report
        .validate_binding(
            &s::selection_revision(
                &catalog,
                &selection,
                "actual-source-role",
                &m,
                &options,
                &batch,
                &choice,
            )
            .unwrap(),
        )
        .unwrap();
    let mut changed = report.clone();
    changed.choice.complement = true;
    assert!(changed.validate_binding(&report.context_revision).is_err());
    changed = report.clone();
    changed.product.product.fragment_index = 1;
    assert!(changed.validate().is_err());
    assert!(report.validate_binding(&"0".repeat(64)).is_err());
}
