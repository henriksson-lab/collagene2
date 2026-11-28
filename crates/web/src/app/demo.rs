use collagene_model::Feature;
use collagene_model::Location;
use collagene_model::Molecule;
use collagene_model::Topology;
use leptos::prelude::*;

pub(super) fn demonstration() -> Molecule {
    let mut seed = 271828u32;
    let mut bases: String = (0..4800)
        .map(|_| {
            seed = seed.wrapping_mul(1664525).wrapping_add(1013904223);
            b"ACGT"[(seed >> 24) as usize % 4] as char
        })
        .collect();
    for (pos, site) in [
        (250, "GAATTC"),
        (1250, "GGATCC"),
        (2700, "AAGCTT"),
        (4050, "GCGGCCGC"),
    ] {
        bases.replace_range(pos..pos + site.len(), site);
    }
    let mut mol = Molecule::dna("pCG_demo", &bases, Topology::Circular).expect("synthetic DNA");
    for (kind, label, start, end, reverse) in [
        ("promoter", "Promoter example", 120, 320, false),
        ("CDS", "Reporter example", 450, 1250, false),
        ("terminator", "Terminator example", 1300, 1420, false),
        ("CDS", "Marker example", 1900, 2750, true),
        (
            "rep_origin",
            "Replication origin example",
            3200,
            3700,
            false,
        ),
    ] {
        let location = Location::range(start, end);
        mol.features.push(Feature {
            kind: kind.into(),
            location: if reverse {
                Location::Complement(Box::new(location))
            } else {
                location
            },
            qualifiers: vec![("label".into(), Some(label.into()))],
            needs_review: false,
        });
    }
    mol.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::Join(vec![Location::range(4650, 4800), Location::range(0, 60)]),
        qualifiers: vec![("label".into(), Some("Across-origin example".into()))],
        needs_review: false,
    });
    mol.metadata
        .insert("collagene:example".into(), "interface-demonstration".into());
    mol.metadata.insert("description".into(),"Synthetic interface demonstration. Annotations illustrate the UI and do not describe a lab-ready construct.".into());
    mol
}
