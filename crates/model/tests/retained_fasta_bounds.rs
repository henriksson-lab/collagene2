//! Retained original bytes have their existing ledger; ordinary text stays strict.
use collagene_model::{Error, Molecule, Topology, bounds};
const MIB: usize = 1024 * 1024;
fn molecule() -> Molecule {
    Molecule::dna("literal byte-ledger fixture", "ACGT", Topology::Linear).unwrap()
}

#[test]
fn exact_known_fasta_key_over_one_mib_is_admitted_and_current_decoders_preserve_it() {
    let mut source = molecule();
    let bytes = "A".repeat(MIB + 4);
    source
        .metadata
        .insert("source:original_fasta_base64".into(), bytes.clone());
    bounds::model_preflight(&[&source]).unwrap();
    let json = serde_json::to_vec(&source).unwrap();
    let decoded: Molecule = serde_json::from_slice(&json).unwrap();
    assert_eq!(decoded, source);
    let wire = postcard::to_allocvec(&source).unwrap();
    let decoded: Molecule = postcard::from_bytes(&wire).unwrap();
    assert_eq!(decoded, source);
    // A similar-looking unknown key is ordinary metadata, not an escape hatch.
    let mut ordinary = molecule();
    ordinary
        .metadata
        .insert("source:original_fasta_base64_extra".into(), bytes);
    assert!(matches!(
        bounds::model_preflight(&[&ordinary]),
        Err(Error::TooLarge)
    ));
    assert!(serde_json::from_slice::<Molecule>(&serde_json::to_vec(&ordinary).unwrap()).is_err());
}

#[test]
fn per_value_original_caps_share_the_text_ledger_across_keys_and_models() {
    let mut first = molecule();
    first
        .metadata
        .insert("source:original_fasta_base64".into(), "A".repeat(3 * MIB));
    let mut second = molecule();
    second.original_genbank = Some("A".repeat(3 * MIB));
    // Six MiB in two valid originals is allowed: four MiB is a per-value cap.
    bounds::model_preflight(&[&first, &second]).unwrap();
    first
        .metadata
        .insert("trace_original_base64".into(), "A".repeat(3 * MIB));
    // Both models remain individually valid, but their nine MiB of originals
    // exceed the shared eight MiB text ledger when admitted together.
    bounds::model_preflight(&[&first]).unwrap();
    bounds::model_preflight(&[&second]).unwrap();
    assert!(matches!(
        bounds::model_preflight(&[&first, &second]),
        Err(Error::TooLarge)
    ));
    let mut single = molecule();
    single.metadata.insert(
        "source:original_fasta_base64".into(),
        "A".repeat(4 * MIB + 4),
    );
    assert!(matches!(
        bounds::model_preflight(&[&single]),
        Err(Error::TooLarge)
    ));
    // Import-scale annotation/native storage retains its separate32MiB profile.
    bounds::annotation_edit_preflight(&[&single]).unwrap();
}

#[test]
fn known_originals_still_consume_the_scientific_eight_mib_shared_text_ledger() {
    let mut source = molecule();
    source
        .metadata
        .insert("source:original_fasta_base64".into(), "A".repeat(3 * MIB));
    for index in 0..4 {
        source
            .metadata
            .insert(format!("ordinary-{index}"), "A".repeat(MIB));
    }
    bounds::model_preflight(&[&source]).unwrap();
    source.metadata.insert("ordinary-4".into(), "A".repeat(MIB));
    assert!(matches!(
        bounds::model_preflight(&[&source]),
        Err(Error::TooLarge)
    ));
}
