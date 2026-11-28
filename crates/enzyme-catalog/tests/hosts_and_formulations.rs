//! Independently authored literal host and formulation facts; no molecule inference.
use collagene_enzyme_catalog::{
    Catalog, Condition, FactValue, Resolution, Sensitivity,
    hosts::{HostPreset, MAX_HOSTS, MethylaseExpectation as M},
    seal,
};
fn known(c: &Catalog, id: &str, condition: &Condition) -> FactValue {
    let Resolution::Known { value, .. } = c.resolve(id, condition, None).unwrap() else {
        panic!("missing literal {id} {condition:?}")
    };
    value.clone()
}
#[test]
fn source_scoped_host_roster_preserves_positive_negative_mixed_and_unknown_expectations() {
    let c = Catalog::builtin().unwrap();
    assert_eq!(c.hosts().len(), 8);
    for (id, code, dam, dcm) in [
        ("neb-5-alpha-c2987", "C2987", M::Active, M::Active),
        ("neb-10-beta-c3019", "C3019", M::Active, M::Active),
        ("neb-stable-c3040", "C3040", M::Active, M::Active),
        ("neb-dam-dcm-c2925", "C2925", M::Deficient, M::Deficient),
        ("neb-bl21-de3-c2527", "C2527", M::Active, M::Deficient),
        ("agilent-jm110-200239", "200239", M::Deficient, M::Deficient),
        (
            "thermo-inv110-c717103",
            "C717103",
            M::Deficient,
            M::Deficient,
        ),
        ("thermo-top10-c404003", "C404003", M::Active, M::Unknown),
    ] {
        let h = c.host(id).unwrap();
        assert_eq!((h.product_code.as_str(), h.dam, h.dcm), (code, dam, dcm));
        assert_eq!(h.cpg, M::Unknown);
        assert_eq!(h.species, "Escherichia coli");
        assert!(!h.vendor.is_empty() && !h.host_scope.is_empty() && !h.citations.is_empty());
        assert!(
            h.qualifiers
                .iter()
                .any(|q| q.contains("no measured plasmid occupancy"))
        );
        c.validate_host(h).unwrap();
    }
    let deficient = c.host("neb-dam-dcm-c2925").unwrap();
    assert!(deficient.host_scope.contains("dcm-6") && deficient.host_scope.contains("dam13::Tn9"));
    assert!(
        deficient
            .qualifiers
            .iter()
            .any(|q| q.contains("M.EcoKI activity"))
    );
    assert!(
        deficient
            .qualifiers
            .iter()
            .any(|q| q.contains("primary cloning/ligation"))
    );
    assert!(c.host("DH5alpha").is_err());
    assert!(c.host("JM110").is_err());
    assert!(c.host("NEB-5-ALPHA-C2987").is_err());
}
#[test]
fn current_host_membership_requires_full_unchanged_facts_not_claimed_identity_or_citation() {
    let c = Catalog::builtin().unwrap();
    let h = c.host("neb-5-alpha-c2987").unwrap();
    let cloned = h.clone();
    c.validate_host(&cloned).unwrap();
    let mut changed = h.clone();
    changed.dam = M::Deficient;
    assert!(c.validate_host(&changed).is_err());
    let mut changed = h.clone();
    changed.cpg = M::Active;
    assert!(c.validate_host(&changed).is_err());
    let mut changed = h.clone();
    changed.citations[0].locator = "changed claimed evidence".into();
    assert!(c.validate_host(&changed).is_err());
    let mut changed = h.clone();
    changed.host_scope = "every laboratory strain".into();
    assert!(c.validate_host(&changed).is_err());
    let mut changed = h.clone();
    changed.qualifiers.clear();
    assert!(changed.validate().is_err());
    let mut changed = h.clone();
    changed.name = "x".repeat(129);
    assert!(c.validate_host(&changed).is_err());
}
#[test]
fn catalogue_hashes_bind_host_assertions_and_require_current_schema_and_hosts_field() {
    let c = Catalog::builtin().unwrap();
    assert_eq!(c.edition().schema_version, 3);
    assert_eq!(c.edition().edition, "2026-11-28.2");
    let mut edited = c.edition().clone();
    edited.hosts[0].dam = M::Active;
    assert!(Catalog::from_edition(edited.clone()).is_err());
    let resealed = Catalog::from_edition(seal(edited).unwrap()).unwrap();
    assert_ne!(resealed.content_hash(), c.content_hash());
    assert!(
        c.validate_host(resealed.host(&c.hosts()[0].id).unwrap())
            .is_err()
    );
    for version in [1, 2, 4] {
        let mut edited = c.edition().clone();
        edited.schema_version = version;
        assert!(seal(edited).is_err());
    }
    let mut missing = serde_json::to_value(c.edition()).unwrap();
    missing.as_object_mut().unwrap().remove("hosts");
    assert!(Catalog::from_json(&serde_json::to_vec(&missing).unwrap()).is_err());
    let round = Catalog::from_json(&c.to_json().unwrap()).unwrap();
    assert_eq!(round.content_hash(), c.content_hash());
}
#[test]
fn host_nested_decode_and_native_proposals_are_bounded_and_postcard_current() {
    let c = Catalog::builtin().unwrap();
    let host = c.hosts()[0].clone();
    let encoded = postcard::to_allocvec(&host).unwrap();
    let decoded: HostPreset = postcard::from_bytes(&encoded).unwrap();
    assert_eq!(decoded, host);
    let mut extra = serde_json::to_value(&host).unwrap();
    extra["fully_methylated"] = true.into();
    assert!(serde_json::from_value::<HostPreset>(extra).is_err());
    let mut oversized = serde_json::to_value(&host).unwrap();
    oversized["qualifiers"] = serde_json::json!(vec!["q"; 17]);
    assert!(serde_json::from_value::<HostPreset>(oversized).is_err());
    let mut oversized = serde_json::to_value(c.edition()).unwrap();
    oversized["hosts"] = serde_json::json!(vec![host.clone(); MAX_HOSTS + 1]);
    assert!(Catalog::from_json(&serde_json::to_vec(&oversized).unwrap()).is_err());
    let mut duplicate = c.edition().clone();
    duplicate.hosts.push(host);
    assert!(seal(duplicate).is_err());
    let mut foreign = c.edition().clone();
    foreign.hosts[0].citations[0].source = "not-in-edition".into();
    assert!(seal(foreign).is_err());
}
#[test]
fn practical_cloning_formulations_keep_product_specific_literal_activity_and_temperature() {
    let c = Catalog::builtin().unwrap();
    assert_eq!(c.edition().enzymes.len(), 40);
    assert_eq!(c.edition().formulations.len(), 8);
    assert_eq!(
        c.edition()
            .formulations
            .iter()
            .map(|f| f.evidence.len())
            .sum::<usize>(),
        64
    );
    for (id, temp, activities) in [
        (
            "R0136",
            37000,
            [
                FactValue::Percent(75),
                FactValue::Percent(100),
                FactValue::Percent(100),
                FactValue::Percent(100),
            ],
        ),
        (
            "R3136",
            37000,
            [
                FactValue::Percent(100),
                FactValue::Percent(50),
                FactValue::Percent(10),
                FactValue::Percent(100),
            ],
        ),
        (
            "R3733",
            37000,
            [
                FactValue::Percent(100),
                FactValue::Percent(100),
                FactValue::Percent(100),
                FactValue::Percent(100),
            ],
        ),
        (
            "R0739",
            55000,
            [
                FactValue::LessThanPercent(10),
                FactValue::Percent(50),
                FactValue::Percent(100),
                FactValue::Percent(25),
            ],
        ),
    ] {
        assert_eq!(
            known(&c, id, &Condition::IncubationTemperature),
            FactValue::TemperatureMilliC(temp)
        );
        for (buffer, value) in [
            "NEBuffer r1.1",
            "NEBuffer r2.1",
            "NEBuffer r3.1",
            "rCutSmart",
        ]
        .into_iter()
        .zip(activities)
        {
            assert_eq!(
                known(
                    &c,
                    id,
                    &Condition::BufferActivity {
                        buffer: buffer.into(),
                        temperature_milli_c: temp
                    }
                ),
                value
            );
        }
    }
    let unknown = Condition::BufferActivity {
        buffer: "rCutSmart".into(),
        temperature_milli_c: 37000,
    };
    assert!(
        matches!(c.resolve("R0739",&unknown,None).unwrap(),Resolution::Unknown(v)if v.is_empty())
    );
    assert_eq!(
        known(
            &c,
            "R3733",
            &Condition::Methylation {
                system: collagene_enzyme_catalog::MethylationSystem::Dcm
            }
        ),
        FactValue::Methylation(Sensitivity::ImpairedSomeOverlappingCombinations)
    );
    assert_eq!(
        known(
            &c,
            "R0739",
            &Condition::Methylation {
                system: collagene_enzyme_catalog::MethylationSystem::CpG
            }
        ),
        FactValue::Methylation(Sensitivity::Blocked)
    );
}
#[test]
fn star_assay_context_and_censoring_remain_documentary_without_universal_predictions() {
    let c = Catalog::builtin().unwrap();
    for (id, amounts) in [
        ("R0136", ["32 units", "4 units"]),
        (
            "R3136",
            [">=4000 units", "Reduced star activity is not absence"],
        ),
    ] {
        let Resolution::Known { evidence, .. } = c
            .resolve(id, &Condition::IncubationTemperature, None)
            .unwrap()
        else {
            unreachable!()
        };
        let qualifiers = evidence
            .iter()
            .flat_map(|e| e.qualifiers.iter())
            .cloned()
            .collect::<Vec<_>>()
            .join(" ");
        assert!(
            qualifiers.contains("lot unspecified")
                && qualifiers.contains("5 percent glycerol")
                && qualifiers.contains("1 microgram/50 microlitres")
                && qualifiers.contains("1 hour at 37 C")
        );
        for value in amounts {
            assert!(qualifiers.contains(value));
        }
        assert!(
            evidence
                .iter()
                .any(|e| e.citations.iter().any(|c| c.source == "neb-bamhi-star-faq"))
        );
    }
    let key = Condition::BufferActivity {
        buffer: "rCutSmart".into(),
        temperature_milli_c: 37000,
    };
    let Resolution::Known { evidence, .. } = c.resolve("R0136", &key, None).unwrap() else {
        unreachable!()
    };
    assert!(
        evidence
            .iter()
            .flat_map(|e| e.qualifiers.iter())
            .any(|q| q.contains("does not establish specificity"))
    );
}
