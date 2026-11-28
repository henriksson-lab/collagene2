//! Reference-document transformations over the canonical annotation projection kernel.
use collagene_model::annotations::{
    MappingSegment, Projection, project_features, project_features_with_topology,
};
use collagene_model::{Alphabet, Error, Feature, Molecule, Sequence, Topology};

const EVIDENCE_GENBANK: &str = "annotation_mapping:original_genbank";
const EVIDENCE_EMBL: &str = "annotation_mapping:original_embl";
fn bad(message: impl Into<String>) -> Error {
    Error::Coordinates(message.into())
}
fn unsupported(message: impl Into<String>) -> Error {
    Error::Unsupported(message.into())
}

fn ensure_document_transform(source: &Molecule) -> Result<(), Error> {
    source.validate()?;
    if source.metadata.contains_key("duplex_model")
        || source
            .metadata
            .get("source:snapgene:chemistry_unmapped")
            .is_some_and(|v| v == "true")
    {
        return Err(unsupported(
            "document carries physical end/nick evidence that cannot be transformed as sequence-only data; reconstruct the originating graph",
        ));
    }
    Ok(())
}
fn derived(
    source: &Molecule,
    sequence: Sequence,
    topology: Topology,
    projection: Projection,
    operation: &str,
) -> Result<Molecule, Error> {
    let mut result = source.clone();
    result.sequence = sequence;
    result.topology = topology;
    result.features = projection.features;
    // Exporters use these fields as interpreted current headers. Move them to explicit source
    // evidence so palindromic/equal-sequence coordinate changes cannot reuse stale reference ranges.
    let source_revision = source.revision();
    let retain = |metadata: &mut std::collections::BTreeMap<String, String>,
                  key: &str,
                  value: String|
     -> Result<(), Error> {
        if let Some(existing) = metadata.get(key) {
            if existing == &value {
                return Ok(());
            }
            let additional = format!("{key}:{source_revision}");
            if metadata.get(&additional).is_some_and(|old| old != &value) {
                return Err(unsupported(
                    "source evidence key collision; no attachment was overwritten",
                ));
            }
            metadata.insert(additional, value);
        } else {
            metadata.insert(key.into(), value);
        }
        Ok(())
    };
    if let Some(original) = result.original_genbank.take() {
        // Imported documents can omit explicit molecule-type metadata. Moving the
        // raw LOCUS into evidence must not erase a known single/mixed-strand restriction.
        if let Some(kind) = original
            .lines()
            .find(|line| line.starts_with("LOCUS"))
            .and_then(|line| {
                line.split_whitespace().find(|token| {
                    let upper = token.to_ascii_uppercase();
                    upper.starts_with("SS-") || upper.starts_with("MS-")
                })
            })
        {
            result
                .metadata
                .insert("genbank.molecule_type".into(), kind.into());
        }
        retain(&mut result.metadata, EVIDENCE_GENBANK, original)?;
    }
    if let Some(original) = result.metadata.remove("embl.original_record") {
        retain(&mut result.metadata, EVIDENCE_EMBL, original)?;
    }
    // The source-bound keyboard edit log addresses its original active header,
    // not this transformed reference. Retain it once as unverified source data;
    // later generations with no new active recipe do not create new archives.
    if let Some(recipe) = result
        .metadata
        .remove(collagene_model::GENBANK_EDIT_RECIPE_KEY)
    {
        if recipe.len() > collagene_model::MAX_GENBANK_EDIT_RECIPE_BYTES {
            return Err(Error::TooLarge);
        }
        retain(
            &mut result.metadata,
            "annotation_mapping:unverified_genbank_edit_recipe",
            recipe,
        )?;
    }
    if let Some(scope) = result.metadata.remove("genbank:citation_scope") {
        retain(
            &mut result.metadata,
            "annotation_mapping:ancestral_genbank_citation_scope",
            scope,
        )?;
    }
    for key in [
        "source_record_index",
        "embl.sequence_version",
        "embl.data_class",
        "embl.division",
    ] {
        result.metadata.remove(key);
    }
    result
        .metadata
        .insert("annotation_mapping:source_revision".into(), source_revision);
    result
        .metadata
        .insert("annotation_mapping:operation".into(), operation.into());
    result
        .metadata
        .insert("annotation_mapping:derived".into(), "true".into());
    if !projection.unmapped_features.is_empty() {
        // Never replace earlier evidence on subsequent transformations.
        let mut evidence: Vec<Feature> =
            if let Some(old) = result.metadata.get("annotation_mapping:unmapped_features") {
                serde_json::from_str(old).map_err(|_| {
                    unsupported("retained unmapped-feature evidence could not be decoded")
                })?
            } else {
                vec![]
            };
        evidence.extend(projection.unmapped_features);
        if evidence.len() > 100_000 {
            return Err(Error::TooLarge);
        }
        result.metadata.insert(
            "annotation_mapping:unmapped_features".into(),
            serde_json::to_string(&evidence).map_err(|e| Error::Format(e.to_string()))?,
        );
    }
    let mut warnings = projection.warnings;
    warnings.push("This derived document retains original headers/attachments as source evidence; they are not remapped current database or chemical assertions.".into());
    result.metadata.insert(
        "annotation_mapping:warnings".into(),
        serde_json::to_string(&warnings).map_err(|e| Error::Format(e.to_string()))?,
    );
    result.validate()?;
    // Retained source evidence must not create a derived document that native save cannot hold.
    crate::protocol::save_document(&result)?;
    Ok(result)
}

pub fn reverse_complement(source: &Molecule) -> Result<Molecule, Error> {
    ensure_document_transform(source)?;
    let n = source.sequence.len();
    let projection = project_features_with_topology(
        source,
        &[MappingSegment {
            source_start: 0,
            source_end: n,
            destination_start: 0,
            reverse: true,
        }],
        n,
        source.topology,
    )?;
    let mut result = derived(
        source,
        source.sequence.reverse_complement()?,
        source.topology,
        projection,
        "reverse complement",
    )?;
    if source.metadata.contains_key("trace_original_base64") {
        // Validate the old view, then update orientation without changing original attachment bytes.
        let trace =
            crate::review::reads::traces::trace_from_molecule(source)?.reverse_complement()?;
        let view = trace.attached_molecule()?;
        for key in [
            "trace_orientation",
            "trace_source_hash",
            "trace_attachment_note",
        ] {
            if let Some(value) = view.metadata.get(key) {
                result.metadata.insert(key.into(), value.clone());
            }
        }
    }
    if source
        .metadata
        .contains_key(crate::formats::oligos::FASTQ_ORIGINAL_KEY)
    {
        crate::formats::oligos::project_fastq_reverse_complement(source, &mut result)?;
    }
    result.validate()?;
    crate::protocol::save_document(&result)?;
    Ok(result)
}

pub fn rotate_origin(source: &Molecule, new_origin: usize) -> Result<Molecule, Error> {
    ensure_document_transform(source)?;
    if source.topology != Topology::Circular {
        return Err(unsupported(
            "setting an origin requires a circular document",
        ));
    }
    let n = source.sequence.len();
    if new_origin >= n {
        return Err(bad("new origin is outside the source sequence"));
    }
    if new_origin == 0 {
        return Ok(source.clone());
    }
    let segments = [
        MappingSegment {
            source_start: new_origin,
            source_end: n,
            destination_start: 0,
            reverse: false,
        },
        MappingSegment {
            source_start: 0,
            source_end: new_origin,
            destination_start: n - new_origin,
            reverse: false,
        },
    ];
    let projection = project_features_with_topology(source, &segments, n, Topology::Circular)?;
    let bases = format!(
        "{}{}",
        &source.sequence.as_str()[new_origin..],
        &source.sequence.as_str()[..new_origin]
    );
    derived(
        source,
        Sequence::new(source.sequence.alphabet(), &bases)?,
        Topology::Circular,
        projection,
        &format!("origin set to source boundary {new_origin}"),
    )
}

/// A half-open interval. Circular end may exceed sequence length by at most one turn;
/// e.g. start=8,end=13 extracts five bases across the origin of a ten-base circle.
pub fn extract_region(source: &Molecule, start: usize, end: usize) -> Result<Molecule, Error> {
    ensure_document_transform(source)?;
    let n = source.sequence.len();
    if start >= n
        || end <= start
        || end - start > n
        || source.topology == Topology::Linear && end > n
    {
        return Err(bad(
            "extraction interval is empty, out of bounds or longer than one turn",
        ));
    }
    let mut segments = vec![MappingSegment {
        source_start: start,
        source_end: end.min(n),
        destination_start: 0,
        reverse: false,
    }];
    if end > n {
        segments.push(MappingSegment {
            source_start: 0,
            source_end: end - n,
            destination_start: n - start,
            reverse: false,
        });
    }
    let projection = project_features(source, &segments, end - start)?;
    let bases = if end <= n {
        source.sequence.as_str()[start..end].to_owned()
    } else {
        format!(
            "{}{}",
            &source.sequence.as_str()[start..],
            &source.sequence.as_str()[..end - n]
        )
    };
    let mut result = derived(
        source,
        Sequence::new(source.sequence.alphabet(), &bases)?,
        Topology::Linear,
        projection,
        &format!("source interval {start}..{end}"),
    )?;
    if source.metadata.contains_key("trace_original_base64") {
        result.metadata.insert("trace_projection".into(),"Original trace is retained as source evidence. A region extraction is not a remapped chromatogram; use the original attachment to inspect its full calls/signals.".into());
    }
    result.name = format!("{} region {}–{}", source.name, start + 1, end);
    result.validate()?;
    crate::protocol::save_document(&result)?;
    Ok(result)
}

/// Change the declared reference topology; this is not a physical ligation/end-repair operation.
/// Linearizing a circle marks features spanning its new free boundary as interrupted.
pub fn change_reference_topology(source: &Molecule, topology: Topology) -> Result<Molecule, Error> {
    source.validate()?;
    if source.topology == topology {
        return Ok(source.clone());
    }
    ensure_document_transform(source)?;
    if topology == Topology::Circular && source.sequence.alphabet() == Alphabet::Protein {
        return Err(unsupported("circular protein documents are not supported"));
    }
    if source.metadata.contains_key("trace_original_base64") {
        return Err(unsupported(
            "trace reference topology cannot be changed; retain the linear trace evidence and create a separate design document",
        ));
    }
    let n = source.sequence.len();
    let projection = project_features_with_topology(
        source,
        &[MappingSegment {
            source_start: 0,
            source_end: n,
            destination_start: 0,
            reverse: false,
        }],
        n,
        topology,
    )?;
    let mut result = derived(
        source,
        source.sequence.clone(),
        topology,
        projection,
        "change declared reference topology",
    )?;
    result.metadata.insert("annotation_mapping:reference_topology_note".into(), "Reference topology declaration only; no physical junction, phosphorylation, nick repair or experimental circularity was established.".into());
    result.validate()?;
    crate::protocol::save_document(&result)?;
    Ok(result)
}

pub mod batch;
pub mod library;
