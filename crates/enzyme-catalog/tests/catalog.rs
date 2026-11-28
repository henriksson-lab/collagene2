use collagene_enzyme_catalog::*;

#[test]
fn smai_chart_temperature_retains_the_extended_digest_exception() {
    let catalog = Catalog::builtin().unwrap();
    let Resolution::Known { value, evidence } = catalog
        .resolve("R0141", &Condition::IncubationTemperature, None)
        .unwrap()
    else {
        panic!("expected the documented formulation-specific recommendation");
    };
    assert_eq!(value, &FactValue::TemperatureMilliC(37_000));
    assert!(evidence.iter().any(|item| {
        item.qualifiers
            .iter()
            .any(|q| q.contains("25 C") && q.contains("longer than one hour"))
    }));
    assert!(
        evidence
            .iter()
            .any(|item| item.citations.iter().any(|citation| {
                catalog.edition().sources.iter().any(|source| {
            source.id == citation.source
                && source.url
                    == "https://www.neb.com/en/faqs/which-incubation-temperature-should-i-use"
        })
            }))
    );
}

fn enzyme<'a>(c: &'a Catalog, n: &str) -> &'a Enzyme {
    match c.lookup_enzyme(n) {
        Lookup::Unique(e) => e,
        _ => panic!("not unique {n}"),
    }
}
fn proposal() -> Edition {
    Catalog::builtin().unwrap().edition().clone()
}
fn citation() -> Vec<Citation> {
    vec![Citation {
        source: "neb-recognition".into(),
        locator: "independent synthetic test fixture; not a curated vendor assertion".into(),
    }]
}
#[test]
fn independently_written_geometry_expectations() {
    let c = Catalog::builtin().unwrap();
    assert_eq!(c.edition().enzymes.len(), 40);
    for (name, motif, top, bottom) in [
        ("EcoRI", "GAATTC", 1, 5),
        ("BamHI", "GGATCC", 1, 5),
        ("HindIII", "AAGCTT", 1, 5),
        ("EcoRV", "GATATC", 3, 3),
        ("SmaI", "CCCGGG", 3, 3),
        ("NotI", "GCGGCCGC", 2, 6),
        ("XbaI", "TCTAGA", 1, 5),
        ("BsaI", "GGTCTC", 7, 11),
        ("PstI", "CTGCAG", 5, 1),
        ("BsmBI", "CGTCTC", 7, 11),
        ("BbsI", "GAAGAC", 8, 12),
        ("SapI", "GCTCTTC", 8, 11),
        ("HinfI", "GANTC", 1, 4),
    ] {
        let g = &enzyme(&c, name).geometry;
        assert_eq!(
            (&*g.recognition, g.top_cut, g.bottom_cut),
            (motif, top, bottom)
        );
    }
}
#[test]
fn biochemical_identities_are_not_aliases() {
    let c = Catalog::builtin().unwrap();
    assert_eq!(
        c.relation(enzyme(&c, "HpaII"), enzyme(&c, "MspI")).unwrap(),
        Relation::SameCutIsoschizomers
    );
    assert_eq!(
        c.relation(enzyme(&c, "SmaI"), enzyme(&c, "XmaI")).unwrap(),
        Relation::Neoschizomers
    );
    assert_eq!(
        c.relation(enzyme(&c, "Acc65I"), enzyme(&c, "KpnI"))
            .unwrap(),
        Relation::Neoschizomers
    );
    assert_eq!(
        c.relation(enzyme(&c, "EcoRI"), enzyme(&c, "BamHI"))
            .unwrap(),
        Relation::DifferentRecognition
    );
    assert!(matches!(c.lookup_enzyme("EcoRI-HF"), Lookup::Missing));
    assert!(matches!(c.lookup_formulation("R3101"), Lookup::Unique(_)));
}
#[test]
fn reverse_descriptor_equivalence() {
    let mut edition = proposal();
    let mut reverse = edition
        .enzymes
        .iter()
        .find(|e| e.name == "BsaI")
        .unwrap()
        .clone();
    reverse.id = "reversefixture".into();
    reverse.name = "Synthetic reverse descriptor".into();
    reverse.geometry = Geometry {
        recognition: "GAGACC".into(),
        top_cut: -5,
        bottom_cut: -1,
    };
    reverse.citations = citation();
    edition.enzymes.push(reverse);
    let c = Catalog::from_edition(seal(edition).unwrap()).unwrap();
    assert_eq!(
        c.relation(enzyme(&c, "BsaI"), enzyme(&c, "reversefixture"))
            .unwrap(),
        Relation::SameCutIsoschizomers
    );
}
#[test]
fn relation_rejects_foreign_and_altered_equal_id_records() {
    let c = Catalog::builtin().unwrap();
    let original = enzyme(&c, "EcoRI");
    let cloned = original.clone();
    assert_eq!(
        c.relation(original, &cloned).unwrap(),
        Relation::SameIdentity
    );
    let mut altered = cloned.clone();
    altered.geometry.top_cut = i32::MIN;
    assert!(c.relation(original, &altered).is_err());
    assert!(c.relation(&altered, &altered).is_err());
    let mut foreign = cloned;
    foreign.id = "foreign".into();
    assert!(c.relation(&foreign, original).is_err());
    let mut next_edition = proposal();
    next_edition.edition = "synthetic-next-edition".into();
    let changed = next_edition
        .enzymes
        .iter_mut()
        .find(|e| e.id == "ecori")
        .unwrap();
    changed.citations[0].locator = "different edition documentary binding".into();
    let next = Catalog::from_edition(seal(next_edition).unwrap()).unwrap();
    assert!(c.relation(original, enzyme(&next, "EcoRI")).is_err());
    assert!(next.relation(enzyme(&next, "EcoRI"), original).is_err());
}
#[test]
fn independent_formulation_percentages() {
    let c = Catalog::builtin().unwrap();
    let key = Condition::BufferActivity {
        buffer: "rCutSmart".into(),
        temperature_milli_c: 37000,
    };
    for (id, percent) in [("R0101", 50), ("R3101", 100)] {
        assert!(
            matches!(c.resolve(id,&key,None).unwrap(),Resolution::Known{value:FactValue::Percent(p),..}if *p==percent)
        );
    }
    assert!(matches!(
        c.resolve("R0141", &Condition::IncubationTemperature, None)
            .unwrap(),
        Resolution::Known {
            value: FactValue::TemperatureMilliC(37000),
            ..
        }
    ));
}
#[test]
fn documentary_methylation_distinctions() {
    let c = Catalog::builtin().unwrap();
    let key = Condition::Methylation {
        system: MethylationSystem::CpG,
    };
    assert!(matches!(
        c.resolve("R0141", &key, None).unwrap(),
        Resolution::Known {
            value: FactValue::Methylation(Sensitivity::Blocked),
            ..
        }
    ));
    assert!(matches!(
        c.resolve("R0180", &key, None).unwrap(),
        Resolution::Known {
            value: FactValue::Methylation(Sensitivity::Impaired),
            ..
        }
    ));
    let missing = Condition::BufferActivity {
        buffer: "unmeasured".into(),
        temperature_milli_c: 37000,
    };
    assert!(
        matches!(c.resolve("R0180",&missing,None).unwrap(),Resolution::Unknown(v)if v.is_empty())
    );
    assert!(c.resolve("HindIII", &key, None).is_err());
}
#[test]
fn conflicts_and_contexts_are_not_silently_collapsed() {
    let mut e = proposal();
    let f = &mut e.formulations[0];
    let id = f.id.clone();
    let mut fact = f
        .evidence
        .iter()
        .find(|v| v.condition == Condition::IncubationTemperature)
        .unwrap()
        .clone();
    fact.value = FactValue::TemperatureMilliC(25000);
    f.evidence.push(fact.clone());
    let c = Catalog::from_edition(seal(e.clone()).unwrap()).unwrap();
    assert!(matches!(
        c.resolve(&id, &Condition::IncubationTemperature, None)
            .unwrap(),
        Resolution::Conflict(_)
    ));
    e.formulations[0].evidence.last_mut().unwrap().context =
        "different historical formulation assay".into();
    let c = Catalog::from_edition(seal(e).unwrap()).unwrap();
    assert!(
        matches!(c.resolve(&id,&Condition::IncubationTemperature,None).unwrap(),Resolution::ContextRequired(v)if v.len()==2)
    );
    assert!(matches!(
        c.resolve(
            &id,
            &Condition::IncubationTemperature,
            Some("vendor recommended incubation")
        )
        .unwrap(),
        Resolution::Known {
            value: FactValue::TemperatureMilliC(37000),
            ..
        }
    ));
}
#[test]
fn explicit_shared_alias_is_ambiguous_in_stable_order() {
    let mut e = proposal();
    for en in &mut e.enzymes[..2] {
        en.aliases.push(Alias {
            name: "synthetic shared name".into(),
            shared: true,
            citations: citation(),
        })
    }
    let c = Catalog::from_edition(seal(e.clone()).unwrap()).unwrap();
    assert!(
        matches!(c.lookup_enzyme(" SYNTHETIC shared NAME "),Lookup::Ambiguous(v)if v.len()==2&&v[0].id<v[1].id)
    );
    e.enzymes[0].aliases[0].shared = false;
    assert!(seal(e).is_err());
}
#[test]
fn duplicate_shadow_and_bad_citation_rejected() {
    let mut e = proposal();
    e.enzymes.push(e.enzymes[0].clone());
    assert!(seal(e).is_err());
    let mut e = proposal();
    e.enzymes[0].aliases.push(Alias {
        name: "EcoRI".into(),
        shared: true,
        citations: citation(),
    });
    assert!(seal(e).is_err());
    let mut e = proposal();
    e.enzymes[0].citations[0].source = "nonexistent".into();
    assert!(seal(e).is_err());
}
#[test]
fn malformed_geometry_never_panics() {
    for top in [i32::MIN, i32::MAX, -1025, 1025] {
        let mut e = proposal();
        e.enzymes[0].geometry.top_cut = top;
        assert!(seal(e).is_err())
    }
    for motif in ["", "GAUTTC", "gaattc", "ACGT!", &"A".repeat(65)] {
        let mut e = proposal();
        e.enzymes[0].geometry.recognition = motif.into();
        assert!(seal(e).is_err())
    }
}
#[test]
fn unknown_does_not_bypass_condition_validation() {
    let mut e = proposal();
    e.formulations[0].evidence[0].condition = Condition::BufferActivity {
        buffer: "".into(),
        temperature_milli_c: i32::MIN,
    };
    e.formulations[0].evidence[0].value = FactValue::Unknown("no data".into());
    assert!(seal(e).is_err());
}
#[test]
fn immutable_checksums_canonical_order_and_tampering() {
    let c = Catalog::builtin().unwrap();
    let mut e = c.edition().clone();
    e.enzymes.reverse();
    e.formulations.reverse();
    e.sources.reverse();
    assert_eq!(
        Catalog::from_edition(e.clone()).unwrap().content_hash(),
        c.content_hash()
    );
    e.enzymes[0].geometry.top_cut += 1;
    assert!(Catalog::from_edition(e).is_err());
    let d = Catalog::from_json(&c.to_json().unwrap()).unwrap();
    assert_eq!(d.content_hash(), c.content_hash());
}
#[test]
fn bounded_strict_json_and_future_schema() {
    for version in [1, 2, 4] {
        let mut e = proposal();
        e.schema_version = version;
        assert!(seal(e).is_err());
    }
    assert!(Catalog::from_json(&vec![b' '; MAX_JSON_BYTES + 1]).is_err());
    let mut json = serde_json::to_value(proposal()).unwrap();
    json["unexpected"] = true.into();
    assert!(Catalog::from_json(&serde_json::to_vec(&json).unwrap()).is_err());
    let mut json = serde_json::to_value(proposal()).unwrap();
    json["enzymes"][0]["geometry"]["unexpected"] = true.into();
    assert!(Catalog::from_json(&serde_json::to_vec(&json).unwrap()).is_err());
}
#[test]
fn search_is_deterministic_bounded_and_case_insensitive() {
    let c = Catalog::builtin().unwrap();
    assert_eq!(enzyme(&c, " eCoRI ").name, "EcoRI");
    let hits = c.search("hpa", 3).unwrap();
    assert_eq!(
        hits.iter().map(|e| e.name.as_str()).collect::<Vec<_>>(),
        ["HpaI", "HpaII"]
    );
    assert!(c.search("x", 1001).is_err());
    assert!(c.search(&"x".repeat(129), 1).is_err());
}
#[test]
fn upstream_hashes_are_optional_not_fabricated() {
    let c = Catalog::builtin().unwrap();
    assert_eq!(
        c.edition()
            .sources
            .iter()
            .filter(|s| s.snapshot.is_some())
            .count(),
        2
    );
    assert!(
        c.edition()
            .sources
            .iter()
            .all(|s| s.facts_blake3.len() == 64)
    );
}

#[test]
fn native_proposals_apply_the_byte_budget_before_validation_or_checksum_work() {
    let mut e = proposal();
    let mut template = e.formulations[0].clone();
    template.aliases.clear();
    let mut fact = template.evidence[0].clone();
    fact.qualifiers = vec!["q".repeat(1024); 16];
    template.evidence = vec![fact; 64];
    // The two legal-size qualifier arrays alone occupy exactly the complete byte
    // budget; JSON syntax and the existing catalog necessarily exceed it. This
    // modest fixture needs about 2 MiB, rather than a giant invalid single string.
    assert_eq!(2 * 64 * 16 * 1024, MAX_JSON_BYTES);
    for index in 0..2 {
        let mut f = template.clone();
        f.id = format!("synthetic-byte-budget-{index}");
        f.name = format!("Synthetic byte budget formulation {index}");
        e.formulations.push(f);
    }
    assert_eq!(
        Catalog::from_edition(e.clone()).unwrap_err().0,
        "edition exceeds serialized byte limit"
    );
    assert_eq!(
        seal(e).unwrap_err().0,
        "edition exceeds serialized byte limit"
    );
}

#[test]
fn oversized_native_alias_names_fail_before_normalization_or_citation_processing() {
    for formulation in [false, true] {
        let mut e = proposal();
        let alias = Alias {
            name: "A".repeat(129),
            shared: false,
            // Deliberately bad citation distinguishes the earlier alias bound.
            citations: vec![Citation {
                source: "missing-source".into(),
                locator: "independent adversarial fixture".into(),
            }],
        };
        if formulation {
            e.formulations[0].aliases.push(alias);
        } else {
            e.enzymes[0].aliases.push(alias);
        }
        assert_eq!(seal(e).unwrap_err().0, "invalid alias name");
    }
}

#[test]
fn source_and_edition_dates_are_real_gregorian_dates_with_century_leap_rules() {
    for valid in [
        "2024-02-29",
        "2000-02-29",
        "2026-04-30",
        "0001-01-01",
        "9999-12-31",
    ] {
        let mut e = proposal();
        e.curated_on = valid.into();
        for source in &mut e.sources {
            source.retrieved_on = valid.into();
        }
        assert!(seal(e).is_ok(), "valid calendar date {valid}");
    }
    for invalid in [
        "2026-02-31",
        "2026-04-31",
        "2100-02-29",
        "1900-02-29",
        "0000-01-01",
        "2026-2-01",
        "2026-01-00",
        "2026-13-01",
        "２０２６-01-01",
    ] {
        let mut e = proposal();
        e.curated_on = invalid.into();
        assert!(seal(e).is_err(), "invalid edition calendar date {invalid}");
        let mut e = proposal();
        e.sources[0].retrieved_on = invalid.into();
        assert!(seal(e).is_err(), "invalid source calendar date {invalid}");
    }
}
