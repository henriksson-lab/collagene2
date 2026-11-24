//! Fresh, explicit VCF/BCF allele-to-DNA checking. File calls are declarations,
//! not proof of build identity, genotype or function.
use crate::{Alphabet, Error, Molecule};
pub use collagene_variants as exchange;
use serde::{Deserialize, Serialize};
const MAX_ALLELE: usize = 65_536;
pub const ENGINE_VERSION: &str = "variant-reference-check-1";
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Options{
    pub source_blake3:String,pub record_index:usize,
    /// Explicit file-contig to actual-DNA mapping, never inferred from a name.
    pub reference_name:String,pub expected_reference_revision:String,
    /// One-based alternative index. Genotype columns never select it implicitly.
    pub alternate_index:usize,
}}
dto! {#[derive(PartialEq,Eq)] pub enum Outcome{
    ReferenceAlleleMatches,ReferenceAlleleMismatch,AmbiguousReference,
    Unsupported{reason:String},
}}
dto! {pub struct Report{
    pub engine_version:String,pub revision:String,pub options:Options,
    pub source_blake3:String,pub actual_reference_revision:String,
    pub reference_name:String,pub start:Option<usize>,pub end:Option<usize>,
    pub declared_reference:String,pub declared_alternative:String,
    pub actual_reference:Option<String>,pub outcome:Outcome,pub interpretation:String,
}}
dto! {pub struct Preview{pub report:Report,pub molecule:Molecule}}
fn bad(s: impl Into<String>) -> Error {
    Error::Unsupported(format!("variant allele check: {}", s.into()))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn concrete(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
}
fn validate(options: &Options, reference: &Molecule) -> Result<(), Error> {
    if !hex(&options.source_blake3)
        || !hex(&options.expected_reference_revision)
        || options.reference_name.is_empty()
        || options.reference_name.len() > 1024
        || options.reference_name.chars().any(char::is_control)
        || options.alternate_index == 0
        || options.alternate_index > 64
    {
        return Err(bad("invalid source/reference/alternative binding"));
    }
    crate::model::bounds::model_preflight(&[reference])?;
    if reference.sequence.alphabet() != Alphabet::Dna {
        return Err(bad("the explicitly mapped reference must be DNA"));
    }
    reference.validate()
}
/// Readonly evidence only. No API accepts a decoded Report as an editing proof.
pub fn check_record(
    collection: &exchange::Collection,
    reference: &Molecule,
    options: &Options,
) -> Result<Report, Error> {
    validate(options, reference)?;
    if collection.source_blake3() != options.source_blake3 {
        return Err(bad("source collection changed"));
    }
    let actual_reference_revision = reference.revision();
    if actual_reference_revision != options.expected_reference_revision {
        return Err(bad("reference revision changed"));
    }
    let record = collection
        .record(options.record_index)
        .ok_or_else(|| bad("record index is unavailable"))?;
    let name = record
        .reference_sequence_name(collection.header())
        .map_err(|e| bad(e.to_string()))?;
    if name != options.reference_name {
        return Err(bad(
            "record contig differs from the explicit reference mapping",
        ));
    }
    let raw_ref = record.reference_bases();
    if raw_ref.len() == 0 || raw_ref.len() > MAX_ALLELE {
        return Err(Error::TooLarge);
    }
    let mut declared_reference = String::with_capacity(raw_ref.len());
    for b in raw_ref.iter().take(MAX_ALLELE + 1) {
        let b = b.map_err(|e| bad(e.to_string()))?;
        if !b.is_ascii() {
            return Err(bad("non-ASCII reference allele"));
        }
        declared_reference.push(b.to_ascii_uppercase() as char);
        if declared_reference.len() > MAX_ALLELE {
            return Err(Error::TooLarge);
        }
    }
    if declared_reference.len() != raw_ref.len() {
        return Err(bad("reference allele length disagrees with its calls"));
    }
    let alternatives = record.alternate_bases();
    let selected = alternatives
        .iter()
        .nth(options.alternate_index - 1)
        .transpose()
        .map_err(|e| bad(e.to_string()))?
        .ok_or_else(|| bad("selected alternative is absent"))?;
    if selected.len() > MAX_ALLELE {
        return Err(Error::TooLarge);
    }
    let declared_alternative = selected.to_ascii_uppercase();
    let start = record
        .variant_start()
        .transpose()
        .map_err(|e| bad(e.to_string()))?
        .map(|p| usize::from(p) - 1);
    let end = start.and_then(|s| s.checked_add(declared_reference.len()));
    let actual_reference = start
        .zip(end)
        .and_then(|(a, b)| reference.sequence.as_str().get(a..b))
        .map(str::to_owned);
    let outcome = if let Some(length) = collection
        .header()
        .contigs()
        .get(name)
        .and_then(|contig| contig.length())
    {
        if length != reference.sequence.len() {
            Outcome::Unsupported{reason:"Declared contig length differs from the actual mapped DNA; no build identity is assumed".into()}
        } else {
            compare(
                &declared_reference,
                &declared_alternative,
                actual_reference.as_deref(),
            )
        }
    } else {
        compare(
            &declared_reference,
            &declared_alternative,
            actual_reference.as_deref(),
        )
    };
    let revision = blake3::hash(
        &serde_json::to_vec(&(ENGINE_VERSION, options, &actual_reference_revision))
            .map_err(|e| bad(e.to_string()))?,
    )
    .to_hex()
    .to_string();
    Ok(Report{engine_version:ENGINE_VERSION.into(),revision,options:options.clone(),source_blake3:collection.source_blake3().into(),actual_reference_revision,reference_name:name.into(),start,end,declared_reference,declared_alternative,actual_reference,outcome,interpretation:"Only the selected concrete REF/ALT and explicitly mapped actual DNA are checked. Names, header identifiers, REF agreement, quality, INFO and genotype columns do not establish genome-build identity, allele presence, ploidy or function. Structural/symbolic/breakend/spanning alleles need separately modeled geometry. No genotype-selected edit is inferred.".into()})
}
fn compare(declared: &str, alternative: &str, actual: Option<&str>) -> Outcome {
    if !concrete(declared) || !concrete(alternative) {
        return Outcome::Unsupported{reason:"Only explicit nonempty A/C/G/T alleles are editable here; symbolic, breakend, spanning, missing and ambiguous alleles remain source evidence".into()};
    }
    let Some(actual) = actual else {
        return Outcome::Unsupported{reason:"Record has no usable forward reference span within this DNA; circular wrapping is not inferred".into()};
    };
    if !concrete(actual) {
        Outcome::AmbiguousReference
    } else if actual == declared {
        Outcome::ReferenceAlleleMatches
    } else {
        Outcome::ReferenceAlleleMismatch
    }
}
/// Recompute against the actual immutable inputs before constructing an edited
/// reference preview. This is design correspondence, not a variant-calling claim.
pub fn preview(
    collection: &exchange::Collection,
    reference: &Molecule,
    options: &Options,
) -> Result<Preview, Error> {
    let report = check_record(collection, reference, options)?;
    if report.outcome != Outcome::ReferenceAlleleMatches {
        return Err(bad(
            "a fully specified matching reference allele is required before creating a preview",
        ));
    }
    let (Some(start), Some(end)) = (report.start, report.end) else {
        return Err(bad("no checked edit span"));
    };
    let recipe = serde_json::to_string(options).map_err(|e| bad(e.to_string()))?;
    if recipe.len() > 4096 {
        return Err(Error::TooLarge);
    }
    let mut molecule = reference.replace(start, end, &report.declared_alternative)?;
    // A selected file allele supplies a design requirement, not prepared DNA.
    // Keep it usable as an independent target without materializing chemistry.
    molecule
        .metadata
        .insert("source:legacy:chemistry_unmapped".into(), "true".into());
    molecule
        .metadata
        .insert("variant:reference_only".into(), "true".into());
    molecule.metadata.insert(
        "variant:unverified_source".into(),
        report.source_blake3.clone(),
    );
    molecule
        .metadata
        .insert("variant:reference_check".into(), report.revision.clone());
    molecule.metadata.insert("variant:interpretation".into(),"Explicit selected-allele design preview; file evidence is not a measured allele or confirmed genotype".into());
    molecule
        .metadata
        .insert("variant:unverified_recipe".into(), recipe);
    molecule.validate()?;
    Ok(Preview { report, molecule })
}

dto! {pub struct ReplayedPreview {
    pub claimed_revision:String,pub literal_dna_and_topology_match:bool,
    pub fresh:Preview,pub interpretation:String,
}}
/// Historical options locate the request, never authorize its result. Both the
/// original private file collection and exact actual reference are required.
/// Names/annotations may have been edited after opening a design; comparison
/// is explicitly limited to literal DNA and topology, not chemical state.
pub fn replay_preview(
    collection: &exchange::Collection,
    reference: &Molecule,
    claimed: &Molecule,
) -> Result<ReplayedPreview, Error> {
    crate::model::bounds::model_preflight(&[reference, claimed])?;
    claimed.validate()?;
    if claimed.sequence.alphabet() != Alphabet::Dna {
        return Err(bad("historical design must be DNA"));
    }
    let recipe=claimed.metadata.get("variant:unverified_recipe").ok_or_else(||bad("historical design has no retained allele recipe; select its original record and alternative explicitly"))?;
    if recipe.len() > 4096 {
        return Err(Error::TooLarge);
    }
    let options: Options = serde_json::from_str(recipe).map_err(|e| bad(e.to_string()))?;
    let fresh = preview(collection, reference, &options)?;
    let matches =
        claimed.topology == fresh.molecule.topology && claimed.sequence == fresh.molecule.sequence;
    Ok(ReplayedPreview{claimed_revision:claimed.revision(),literal_dna_and_topology_match:matches,fresh,interpretation:"Fresh source-file and actual-reference reconstruction only. Historical recipe text and a matching designed sequence do not restore genotype, biological function, prepared chemistry or earlier verification authority.".into()})
}
