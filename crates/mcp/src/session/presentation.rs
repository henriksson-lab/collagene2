//! Bounded summaries of immutable session data; no new scientific authority.
use collagene_core::formats::alignment as alignment_formats;
use collagene_model::Molecule;
use serde_json::{Value, json};
use std::collections::BTreeMap;

pub(crate) fn summary(id: &str, m: &Molecule) -> Value {
    let source = m
        .metadata
        .iter()
        .filter(|(k, _)| {
            [
                "source:provider",
                "source:canonical_url",
                "source:identifier",
                "source:accession_version",
                "source:retrieved_at",
                "source:content_blake3",
                "source:format",
                "source:sequence_id",
                "source:addgene_origin",
                "source:addgene_completeness",
            ]
            .contains(&k.as_str())
        })
        .collect::<Vec<_>>();
    json!({"id":id,"revision":id,"name":m.name.chars().take(256).collect::<String>(),"length":m.sequence.len(),"alphabet":m.sequence.alphabet(),"topology":m.topology,"feature_count":m.features.len(),"source":source.iter().map(|(k,v)|(k.as_str(),v.chars().take(256).collect::<String>())).collect::<BTreeMap<_,_>>(),"source_values_truncated":source.iter().any(|(_,v)|v.chars().count()>256)})
}

pub(crate) fn alignment_summary(
    id: &str,
    document: &alignment_formats::AlignmentDocument,
) -> Value {
    json!({"alignment":{"id":id,"revision":id},"name":document.name.chars().take(256).collect::<String>(),"alphabet":document.alphabet,"row_count":document.rows.len(),"column_count":document.column_count(),"source_format":document.source.format,"source_blake3":document.source.blake3,"rows":document.rows.iter().take(100).map(|r|json!({"id":r.id,"identifier":r.identifier.chars().take(256).collect::<String>(),"description":r.description.chars().take(256).collect::<String>(),"text_preview_truncated":r.identifier.chars().count()>256||r.description.chars().count()>256})).collect::<Vec<_>>(),"rows_preview_truncated":document.rows.len()>100,"persistence":"session only","semantics":"Imported alignment evidence, not a newly computed alignment or verified external reference coordinate map."})
}
