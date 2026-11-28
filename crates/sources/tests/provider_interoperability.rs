//! Captured public reference records exercise the actual import pipeline offline.
//! The separate live client check is ignored until application HTTP acceptance.
use collagene_core::{
    formats::sequence as formats,
    protocol,
    sequence::analysis::{Enzyme, digest_lengths, restriction_sites},
};
use collagene_model::{Alphabet, Topology};
use collagene_sources::{Config, SourceClient, import_download, parse_source};
use serde_json::{Value, json};

fn validate(report: &formats::ImportReport, provider: &str) -> Result<(String, Value), String> {
    let [molecule] = report.molecules.as_slice() else {
        return Err("complete accession import did not return exactly one record".into());
    };
    if molecule.sequence.alphabet() != Alphabet::Dna
        || molecule.sequence.len() != 5386
        || molecule.topology != Topology::Circular
        || molecule.features.is_empty()
    {
        return Err("J02482.1 identity/topology/annotation evidence differs".into());
    }
    for (key, expected) in [
        ("source:provider", provider),
        ("source:accession_version", "J02482.1"),
    ] {
        if molecule.metadata.get(key).map(String::as_str) != Some(expected) {
            return Err(format!(
                "returned {key} differs from the requested provider/record"
            ));
        }
    }
    let native = protocol::save_document(molecule).map_err(|e| e.to_string())?;
    if protocol::load_document(&native).map_err(|e| e.to_string())? != *molecule {
        return Err("native save/reopen changed the imported document".into());
    }
    let genbank = formats::export_genbank(molecule).map_err(|e| e.to_string())?;
    let reopened = formats::import_genbank(&genbank).map_err(|e| e.to_string())?;
    if reopened.molecules.len() != 1
        || reopened.molecules[0].sequence != molecule.sequence
        || reopened.molecules[0].topology != molecule.topology
        || reopened.molecules[0].features != molecule.features
    {
        return Err("annotated GenBank export changed the sequence/topology/features".into());
    }
    // Independent reference lengths: NEB N3026 product specification, not values
    // produced by the restriction implementation under test. This checks the
    // reference sequence's hypothetical complete dsDNA cut geometry, not physical
    // material state of the phage record or laboratory cleavage success.
    // https://www.neb.com/en-us/products/n3026-x174-dna-haeiii-digest
    let enzyme = Enzyme {
        name: "HaeIII".into(),
        recognition: "GGCC".into(),
        top_cut: 2,
        bottom_cut: 2,
    };
    let sites = restriction_sites(molecule, &enzyme).map_err(|e| e.to_string())?;
    let mut lengths = digest_lengths(molecule, &sites).map_err(|e| e.to_string())?;
    lengths.sort_unstable();
    if lengths != [72, 118, 194, 234, 271, 281, 310, 603, 872, 1078, 1353] {
        return Err(format!("reference HaeIII lengths differ: {lengths:?}"));
    }
    let sequence_hash = blake3::hash(molecule.sequence.as_str().as_bytes())
        .to_hex()
        .to_string();
    Ok((
        sequence_hash.clone(),
        json!({
            "status":"passed", "provider":provider, "accession_version":"J02482.1",
            "canonical_url":molecule.metadata.get("source:canonical_url"),
            "retrieved_at":molecule.metadata.get("source:retrieved_at"),
            "content_blake3":molecule.metadata.get("source:content_blake3"),
            "sequence_blake3":sequence_hash, "length":molecule.sequence.len(),
            "topology":"circular", "features":molecule.features.len(),
            "native_roundtrip_equal":true, "genbank_sequence_topology_features_roundtrip_equal":true,
            "haeiii_reference_lengths":lengths, "import_warnings":report.warnings,
            "physical_material_authority":false,
        }),
    ))
}

#[test]
fn captured_public_records_preserve_annotations_and_independent_reference_digest() {
    let mut first = None;
    for (input, provider, bytes, retrieved_at) in [
        (
            "J02482.1",
            "NCBI nuccore",
            include_bytes!("data/public_records/J02482.1.ncbi.gb").as_slice(),
            "2026-11-28T08:06:13Z",
        ),
        (
            "ddbj:J02482.1",
            "DDBJ",
            include_bytes!("data/public_records/J02482.1.ddbj.gb").as_slice(),
            "2026-11-28T08:06:35Z",
        ),
    ] {
        let report = import_download(&parse_source(input).unwrap(), bytes, retrieved_at).unwrap();
        let (sequence, evidence) = validate(&report, provider).unwrap();
        if let Some(first) = &first {
            assert_eq!(first, &sequence);
        } else {
            first = Some(sequence);
        }
        println!("CAPTURED_SOURCE_EVIDENCE {evidence}");
    }
}

#[tokio::test]
#[ignore = "manual public NCBI/ENA/DDBJ interoperability; provider availability is not a routine CI gate"]
async fn complete_public_records_preserve_identity_roundtrips_and_independent_digest_geometry() {
    // Deliberately anonymous. Config::from_env would consult unrelated user keys.
    let client = SourceClient::new(Config::default()).unwrap();
    let mut expected_sequence = None;
    let mut failures = vec![];
    for (source, provider) in [
        ("J02482.1", "NCBI nuccore"),
        ("ena:J02482.1", "ENA"),
        ("ddbj:J02482.1", "DDBJ"),
    ] {
        let result = match client.import_url(source).await {
            Ok(report) => validate(&report, provider),
            Err(error) => Err(error.to_string()),
        };
        match result {
            Ok((sequence, mut evidence)) => {
                if expected_sequence
                    .as_ref()
                    .is_some_and(|first| first != &sequence)
                {
                    failures.push(format!(
                        "{provider}: normalized DNA differs across providers"
                    ));
                    evidence["status"] = "failed".into();
                    evidence["error"] = "normalized DNA differs across providers".into();
                } else {
                    expected_sequence = Some(sequence);
                }
                println!("LIVE_SOURCE_EVIDENCE {evidence}");
            }
            Err(error) => {
                println!(
                    "LIVE_SOURCE_EVIDENCE {}",
                    json!({"status":"failed","provider":provider,"error":error})
                );
                failures.push(format!("{provider}: {error}"));
            }
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
