//! Annotation-preserving reference transfers. These operations declare sequence
//! intent, not molecular cleavage, ligation or prepared physical material.
use crate::sequence::annotations::{self, MappingSegment};
use collagene_model::{Alphabet, Error, Feature, Location, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::io::Write;

pub const ENGINE_VERSION: &str = "reference-transfer-2";
pub const MAX_PARTS: usize = 4096;
const MAX_BYTES: usize = 32 << 20;
const EVIDENCE_CHUNK: usize = 512 << 10;
const REFERENCE_MODEL: &str =
    "Reference-only annotated copy/paste; physical material state unavailable";
macro_rules! dto {($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Part {
    pub start: usize,
    pub end: usize,
    pub reverse: bool,
}
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Selection {
    pub parts: Vec<Part>,
    pub reverse: bool,
}
dto! {pub struct RequestCopy {pub source:Molecule,pub selection:Selection}}
dto! {pub struct Fragment {
    pub context_revision:String,pub source_revision:String,pub source_length:usize,
    pub source_alphabet:Alphabet,pub selection:Selection,pub molecule:Molecule,
    pub warnings:Vec<String>,pub unmapped_features:Vec<Feature>,pub dropped_features:usize,
}}
dto! {pub struct RequestPaste {
    pub source:Molecule,pub destination:Molecule,pub selection:Selection,
    /// Expected copy_revision, not an imported fragment's asserted sequence hash.
    pub fragment_revision:String,pub start:usize,pub end:usize,
}}
dto! {pub struct Report {
    pub context_revision:String,pub source_revision:String,pub destination_revision:String,
    pub fragment_revision:String,pub selection:Selection,pub start:usize,pub end:usize,
    pub inserted_start:usize,pub inserted_end:usize,pub product:Molecule,
    pub copied_features:usize,pub unmapped_features:Vec<Feature>,pub dropped_features:usize,pub warnings:Vec<String>,
}}
fn bad(value: &str) -> Error {
    Error::Unsupported(format!("reference transfer: {value}"))
}
fn digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
struct Bytes(Vec<u8>);
impl Write for Bytes {
    fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
        if self
            .0
            .len()
            .checked_add(bytes.len())
            .is_none_or(|n| n > MAX_BYTES)
        {
            return Err(std::io::Error::other("reference evidence bound"));
        }
        self.0.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn encoded(value: &impl Serialize) -> Result<String, Error> {
    let mut output = Bytes(Vec::new());
    serde_json::to_writer(&mut output, value).map_err(|_| Error::TooLarge)?;
    String::from_utf8(output.0).map_err(|_| bad("invalid evidence encoding"))
}
fn retain(molecule: &mut Molecule, key: String, value: String) -> Result<(), Error> {
    let mut insert = |key: String, value: String| -> Result<(), Error> {
        if molecule.metadata.get(&key).is_some_and(|old| old != &value) {
            return Err(bad("source evidence identity collision"));
        }
        molecule.metadata.insert(key, value);
        Ok(())
    };
    if value.len() <= EVIDENCE_CHUNK {
        return insert(key, value);
    }
    let mut offset = 0;
    let mut count = 0;
    while offset < value.len() {
        let mut end = (offset + EVIDENCE_CHUNK).min(value.len());
        while !value.is_char_boundary(end) {
            end -= 1;
        }
        insert(format!("{key}:part-{count}"), value[offset..end].into())?;
        offset = end;
        count += 1;
    }
    insert(
        key,
        format!("UTF-8 evidence retained in {count} ordered :part-N entries"),
    )
}
fn preflight(source: &Molecule) -> Result<(), Error> {
    collagene_model::bounds::model_preflight(&[source])?;
    source.validate()?;
    crate::encoding::bounded_json::hash(source)?;
    Ok(())
}
impl Selection {
    /// Same exact biological extraction grammar as Location::extract. A complement
    /// reverses child order as well as complementing each ordered part.
    pub fn from_location(location: &Location, source_length: usize) -> Result<Self, Error> {
        fn visit(
            location: &Location,
            reverse: bool,
            depth: usize,
            parts: &mut Vec<Part>,
        ) -> Result<(), Error> {
            if depth > collagene_model::bounds::MAX_DEPTH || parts.len() >= MAX_PARTS {
                return Err(Error::TooLarge);
            }
            match location {
                Location::Range {
                    start,
                    end,
                    before: false,
                    after: false,
                } => parts.push(Part {
                    start: *start,
                    end: *end,
                    reverse,
                }),
                Location::Complement(inner) => visit(inner, !reverse, depth + 1, parts)?,
                Location::Join(children) => {
                    if children.len() > MAX_PARTS {
                        return Err(Error::TooLarge);
                    }
                    if reverse {
                        for child in children.iter().rev() {
                            visit(child, reverse, depth + 1, parts)?;
                        }
                    } else {
                        for child in children {
                            visit(child, reverse, depth + 1, parts)?;
                        }
                    }
                }
                _ => {
                    return Err(bad(
                        "selection requires exact ranges/joins/complements; fuzzy, ordered, between-base and unresolved locations are not flattened",
                    ));
                }
            }
            Ok(())
        }
        let mut parts = Vec::new();
        visit(location, false, 0, &mut parts)?;
        let selection = Self {
            parts,
            reverse: false,
        };
        selection.validate(source_length, Alphabet::Dna)?;
        location.validate(source_length)?;
        Ok(selection)
    }
    fn validate(&self, length: usize, alphabet: Alphabet) -> Result<usize, Error> {
        if self.parts.is_empty() || self.parts.len() > MAX_PARTS {
            return Err(bad("select 1..4096 nonempty source parts"));
        }
        if alphabet == Alphabet::Protein && (self.reverse || self.parts.iter().any(|p| p.reverse)) {
            return Err(bad("proteins have no reverse complement"));
        }
        let mut total = 0usize;
        for part in &self.parts {
            if part.start >= part.end || part.end > length {
                return Err(bad("copy interval is empty or outside its source"));
            }
            total = total
                .checked_add(part.end - part.start)
                .filter(|n| *n <= collagene_model::MAX_SEQUENCE_LEN)
                .ok_or(Error::TooLarge)?;
        }
        let mut sorted = self
            .parts
            .iter()
            .map(|p| (p.start, p.end))
            .collect::<Vec<_>>();
        sorted.sort_unstable();
        if sorted.windows(2).any(|pair| pair[0].1 > pair[1].0) {
            return Err(bad(
                "repeated or overlapping source parts require a separately declared copy; no one-to-one annotation map was invented",
            ));
        }
        Ok(total)
    }
    fn maps(&self, source: &Molecule) -> Result<(Vec<MappingSegment>, String), Error> {
        let length = self.validate(source.sequence.len(), source.sequence.alphabet())?;
        let mut maps = Vec::with_capacity(self.parts.len());
        let mut letters = String::with_capacity(length);
        let effective = if self.reverse {
            self.parts.iter().rev().collect::<Vec<_>>()
        } else {
            self.parts.iter().collect()
        };
        for part in effective {
            let reverse = part.reverse ^ self.reverse;
            let bases = &source.sequence.as_str()[part.start..part.end];
            maps.push(MappingSegment {
                source_start: part.start,
                source_end: part.end,
                destination_start: letters.len(),
                reverse,
            });
            if reverse {
                letters.push_str(
                    Sequence::new(source.sequence.alphabet(), bases)?
                        .reverse_complement()?
                        .as_str(),
                );
            } else {
                letters.push_str(bases);
            }
        }
        Ok((maps, letters))
    }
}
pub fn copy_revision(request: &RequestCopy) -> Result<String, Error> {
    preflight(&request.source)?;
    request.selection.validate(
        request.source.sequence.len(),
        request.source.sequence.alphabet(),
    )?;
    crate::encoding::bounded_json::hash(&(
        ENGINE_VERSION,
        "copy",
        &request.source,
        &request.selection,
    ))
}
fn mark_reference(molecule: &mut Molecule, role: &str, source_revision: &str) -> Result<(), Error> {
    let active = molecule
        .metadata
        .iter()
        .filter(|(key, _)| {
            key.as_str() == "duplex_model" || key.starts_with("trace_") || key.starts_with("fastq:")
        })
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect::<std::collections::BTreeMap<_, _>>();
    if !active.is_empty() {
        retain(
            molecule,
            format!("reference_transfer:{role}_attachments:{source_revision}"),
            encoded(&active)?,
        )?;
        for key in active.keys() {
            molecule.metadata.remove(key);
        }
    }
    molecule
        .metadata
        .insert("duplex_model".into(), REFERENCE_MODEL.into());
    molecule
        .metadata
        .insert("reference_transfer:reference_only".into(), "true".into());
    molecule
        .metadata
        .insert("reference_transfer:engine".into(), ENGINE_VERSION.into());
    Ok(())
}
fn preserve_source(molecule: &mut Molecule, source: &Molecule, role: &str) -> Result<(), Error> {
    let revision = source.revision();
    if !source.features.is_empty() {
        retain(
            molecule,
            format!("reference_transfer:{role}_features:{revision}"),
            encoded(&source.features)?,
        )?;
    }
    if let Some(original) = &source.original_genbank {
        retain(
            molecule,
            format!("reference_transfer:{role}_genbank:{revision}"),
            original.clone(),
        )?;
    }
    // Raw chemical/trace descriptors remain ancestral evidence only. No copied
    // descriptor is installed as current prepared chemistry or remapped calls.
    let evidence = source
        .metadata
        .iter()
        .filter(|(key, _)| {
            key.as_str() == "duplex_model"
                || key.as_str() == "source:snapgene:chemistry_unmapped"
                || key.starts_with("trace_")
                || key.starts_with("fastq:")
                || key.as_str() == "embl.original_record"
        })
        .map(|(k, v)| (k.clone(), v.clone()))
        .collect::<std::collections::BTreeMap<_, _>>();
    if !evidence.is_empty() {
        retain(
            molecule,
            format!("reference_transfer:{role}_metadata:{revision}"),
            encoded(&evidence)?,
        )?;
    }
    Ok(())
}
fn source_citations(features: &mut [Feature]) -> usize {
    let mut changed = 0;
    for feature in features {
        let mut found = false;
        for (key, _) in &mut feature.qualifiers {
            if key.eq_ignore_ascii_case("citation") {
                *key = "source_citation".into();
                found = true;
            }
        }
        if found {
            feature.needs_review = true;
            changed += 1;
            feature.qualifiers.push(("note".into(),Some("Source citation numbers are retained as source_citation evidence; they do not address destination references.".into())));
        }
    }
    changed
}
pub fn copy(request: &RequestCopy) -> Result<Fragment, Error> {
    let context_revision = copy_revision(request)?;
    let source = &request.source;
    let source_revision = source.revision();
    let (maps, letters) = request.selection.maps(source)?;
    let projection = annotations::project_features(source, &maps, letters.len())?;
    let mut molecule = Molecule {
        name: source.name.clone(),
        sequence: Sequence::new(source.sequence.alphabet(), &letters)?,
        topology: Topology::Linear,
        features: projection.features,
        metadata: Default::default(),
        original_genbank: None,
    };
    let mut warnings = projection.warnings;
    if source_citations(&mut molecule.features) > 0 {
        warnings.push("Copied source citation numbers remain source evidence, not destination reference numbers.".into());
    }
    preserve_source(&mut molecule, source, "source")?;
    mark_reference(&mut molecule, "source", &source_revision)?;
    molecule.metadata.insert(
        "reference_transfer:source_revision".into(),
        source_revision.clone(),
    );
    molecule.metadata.insert(
        "reference_transfer:copy_context".into(),
        context_revision.clone(),
    );
    molecule.metadata.insert(
        "reference_transfer:selection".into(),
        encoded(&request.selection)?,
    );
    warnings.push("Annotated reference copy only; no cleavage, end chemistry, synthesis, repair, ligation or experimental success is established.".into());
    let fragment = Fragment {
        context_revision,
        source_revision,
        source_length: source.sequence.len(),
        source_alphabet: source.sequence.alphabet(),
        selection: request.selection.clone(),
        molecule,
        warnings,
        unmapped_features: projection.unmapped_features,
        dropped_features: projection.dropped_features,
    };
    fragment.validate_binding(&fragment.context_revision)?;
    Ok(fragment)
}
pub fn paste_revision(request: &RequestPaste) -> Result<String, Error> {
    collagene_model::bounds::model_preflight(&[&request.source, &request.destination])?;
    request.source.validate()?;
    request.destination.validate()?;
    request.selection.validate(
        request.source.sequence.len(),
        request.source.sequence.alphabet(),
    )?;
    if request.destination.features.iter().any(|feature| {
        feature.qualifiers.iter().any(|(key, _)| {
            ["transl_except", "anticodon", "rpt_unit_range", "coded_by"]
                .iter()
                .any(|known| key.eq_ignore_ascii_case(known))
        })
    }) {
        return Err(bad(
            "destination coordinate-bearing qualifiers cannot yet be safely remapped",
        ));
    }
    if request.source.sequence.alphabet() != request.destination.sequence.alphabet() {
        return Err(bad("source and destination require matching alphabets"));
    }
    if request.start > request.end
        || request.end > request.destination.sequence.len()
        || !digest(&request.fragment_revision)
    {
        return Err(bad(
            "paste requires an explicit in-range destination interval/caret and copied context",
        ));
    }
    crate::encoding::bounded_json::hash(&(ENGINE_VERSION, "paste", request))
}
pub fn paste(request: &RequestPaste) -> Result<Report, Error> {
    let context_revision = paste_revision(request)?;
    let fragment = copy(&RequestCopy {
        source: request.source.clone(),
        selection: request.selection.clone(),
    })?;
    if fragment.context_revision != request.fragment_revision {
        return Err(bad(
            "copied source or declaration changed; copy again before pasting",
        ));
    }
    let destination_revision = request.destination.revision();
    let mut product = request.destination.replace(
        request.start,
        request.end,
        fragment.molecule.sequence.as_str(),
    )?;
    // Preserve native edit geometry while handling case variants of imported
    // translation assertions consistently with the annotation projector.
    if product.features.len() != request.destination.features.len() {
        return Err(bad(
            "native reference edit changed feature correspondence unexpectedly",
        ));
    }
    for (original, edited) in request
        .destination
        .features
        .iter()
        .zip(&mut product.features)
    {
        let touched = original.location.spans().iter().any(|&(start, end, _)| {
            if request.start == request.end {
                start < request.start && request.start < end
            } else {
                start < request.end && request.start < end
            }
        });
        if touched {
            edited.qualifiers.retain(|(key, _)| {
                !key.eq_ignore_ascii_case("translation") && !key.eq_ignore_ascii_case("codon_start")
            });
        }
    }
    // Project from the original source into the final topology. The clipboard
    // preview is linear and can mark an origin-spanning circular feature broken;
    // carrying that intermediate annotation forward would lose valid continuity
    // when the declared destination is circular again.
    let (mut maps, _) = request.selection.maps(&request.source)?;
    for segment in &mut maps {
        segment.destination_start += request.start;
    }
    let projection = annotations::project_features_with_topology(
        &request.source,
        &maps,
        product.sequence.len(),
        product.topology,
    )?;
    let mut copied = projection.features;
    let source_citation_count = source_citations(&mut copied);
    let copied_features = copied.len();
    product.features.extend(copied);
    preserve_source(&mut product, &request.source, "source")?;
    preserve_source(&mut product, &request.destination, "destination")?;
    mark_reference(&mut product, "destination", &destination_revision)?;
    product.metadata.insert(
        "reference_transfer:source_revision".into(),
        fragment.source_revision.clone(),
    );
    product.metadata.insert(
        "reference_transfer:destination_revision".into(),
        destination_revision.clone(),
    );
    product.metadata.insert(
        "reference_transfer:copy_context".into(),
        fragment.context_revision.clone(),
    );
    product.metadata.insert(
        "reference_transfer:paste_context".into(),
        context_revision.clone(),
    );
    product.metadata.insert(
        "reference_transfer:selection".into(),
        encoded(&request.selection)?,
    );
    let unmapped_features = projection.unmapped_features;
    if !unmapped_features.is_empty() {
        retain(
            &mut product,
            format!(
                "reference_transfer:unmapped_features:{}",
                fragment.source_revision
            ),
            encoded(&unmapped_features)?,
        )?;
    }
    let mut warnings = projection.warnings;
    if source_citation_count > 0 {
        warnings.push("Copied source citation numbers remain source evidence, not destination reference numbers.".into());
    }
    warnings.push("Reference paste only; no cleavage, end chemistry, synthesis, repair, ligation or experimental success is established.".into());
    warnings.push("Destination reference annotations use the native edit transform; touched annotations require review and stale translation assertions are removed. Original source/destination features remain evidence.".into());
    let report = Report {
        context_revision,
        source_revision: fragment.source_revision,
        destination_revision,
        fragment_revision: fragment.context_revision,
        selection: request.selection.clone(),
        start: request.start,
        end: request.end,
        inserted_start: request.start,
        inserted_end: request.start + fragment.molecule.sequence.len(),
        product,
        copied_features,
        unmapped_features,
        dropped_features: projection.dropped_features,
        warnings,
    };
    report.validate_binding(&report.context_revision)?;
    Ok(report)
}
impl Fragment {
    /// Readonly structure guard; paste reconstructs from the actual original source.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        collagene_model::bounds::model_preflight(&[&self.molecule])?;
        collagene_model::bounds::feature_preflight(&[&self.unmapped_features])?;
        if self.source_length > collagene_model::MAX_SEQUENCE_LEN {
            return Err(Error::TooLarge);
        }
        let length = self
            .selection
            .validate(self.source_length, self.source_alphabet)?;
        self.molecule.validate()?;
        if !digest(expected)
            || self.context_revision != expected
            || !digest(&self.source_revision)
            || self.molecule.sequence.len() != length
            || self.molecule.sequence.alphabet() != self.source_alphabet
            || self.molecule.topology != Topology::Linear
            || self
                .molecule
                .metadata
                .get("duplex_model")
                .map(String::as_str)
                != Some(REFERENCE_MODEL)
            || self
                .molecule
                .metadata
                .get("reference_transfer:copy_context")
                != Some(&self.context_revision)
        {
            return Err(bad("copied reference receipt differs from its declaration"));
        }
        crate::encoding::bounded_json::hash(self)?;
        Ok(())
    }
}
impl Report {
    /// Readonly structure guard, not authority for imported clipboard content.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        collagene_model::bounds::model_preflight(&[&self.product])?;
        collagene_model::bounds::feature_preflight(&[&self.unmapped_features])?;
        self.product.validate()?;
        let source_bound = self
            .selection
            .parts
            .iter()
            .map(|part| part.end)
            .max()
            .unwrap_or(0);
        if source_bound > collagene_model::MAX_SEQUENCE_LEN {
            return Err(Error::TooLarge);
        }
        let inserted_length = self
            .selection
            .validate(source_bound, self.product.sequence.alphabet())?;
        if self.inserted_end.checked_sub(self.inserted_start) != Some(inserted_length) {
            return Err(bad("inserted length differs from the copied declaration"));
        }
        if !digest(expected)
            || self.context_revision != expected
            || !digest(&self.source_revision)
            || !digest(&self.destination_revision)
            || !digest(&self.fragment_revision)
            || self.start > self.end
            || self.inserted_start != self.start
            || self.inserted_start >= self.inserted_end
            || self.inserted_end > self.product.sequence.len()
            || self.copied_features > self.product.features.len()
            || self
                .product
                .metadata
                .get("duplex_model")
                .map(String::as_str)
                != Some(REFERENCE_MODEL)
            || self
                .product
                .metadata
                .get("reference_transfer:paste_context")
                != Some(&self.context_revision)
            || self.product.metadata.get("reference_transfer:copy_context")
                != Some(&self.fragment_revision)
        {
            return Err(bad("pasted reference receipt differs from its declaration"));
        }
        crate::encoding::bounded_json::hash(self)?;
        Ok(())
    }
}
