use crate::server::response::bad;
use crate::session::genome_store;
use rmcp::ErrorData;
use serde_json::Value;
use serde_json::json;
use std::collections::BTreeMap;

pub(crate) struct NativeSearch {
    pub(crate) genome: collagene_core::tracks::genome_session::GenomeRef,
    pub(crate) state: collagene_core::tracks::crispr::GenomeSearchSession,
}
#[derive(Default)]
pub(crate) struct GenomeResources {
    pub(crate) native: Option<genome_store::NativeGenomes>,
    pub(crate) searches: BTreeMap<String, NativeSearch>,
    pub(crate) serial: u64,
}
impl GenomeResources {
    pub(crate) fn native(&mut self) -> Result<&mut genome_store::NativeGenomes, ErrorData> {
        if self.native.is_none() {
            self.native = Some(genome_store::NativeGenomes::new().map_err(bad)?);
        }
        Ok(self.native.as_mut().unwrap())
    }
    pub(crate) fn existing(&mut self) -> Result<&mut genome_store::NativeGenomes, ErrorData> {
        self.native
            .as_mut()
            .ok_or_else(|| bad("no live large-reference session"))
    }
}
pub(crate) fn native_ref(
    d: &genome_store::NativeDataset,
) -> collagene_core::tracks::genome_session::GenomeRef {
    collagene_core::tracks::genome_session::GenomeRef {
        id: d.id.clone(),
        manifest_hash: d.reference.manifest_hash.clone(),
        content_revision: d.reference.manifest.content_revision.clone(),
    }
}
pub(crate) fn native_summary(d: &genome_store::NativeDataset) -> Value {
    let m = &d.reference.manifest;
    json!({"genome":native_ref(d),"name":m.build.name.chars().take(256).collect::<String>(),"assembly_accession":m.build.assembly_accession,"version":m.build.version,"contig_count":m.contigs.len(),"total_bases":m.total_bases,"input_bytes":d.input_bytes,"file_length":d.reference.file_length,"semantic_verified_at_import":d.semantic_verified,"expected_components_complete":m.expected_complete,"missing_expected_count":m.missing_expected.len(),"raw_source_hash":m.raw_source_hash,"storage":"Session-owned managed scratchfile; no arbitrary filesystem path or persistent workspace write."})
}
