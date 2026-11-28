//! Target-free native PCR, reciprocal extension, catalog-digest and end-repair derivations. Opaque, bounded recipe
//! packages break typed recursion; only fresh native execution yields material.
use crate::documents::history::transition::{AnalysisOrigin, AnalysisProducer, Transition};
pub mod digest_export;
pub mod reaction_conditions;
use crate::{
    chemistry::{
        self,
        tracked_duplex::{self as modified_duplex, ModificationSummary, TrackedDuplex},
    },
    cloning::restriction::{
        catalog::{self as catalog_digest, CatalogSelection},
        inspection::{self as restriction_inspection, InspectionOptions},
    },
    documents::{
        history::DocumentHistory,
        workspace::{Workspace, WorkspaceAttachment},
    },
    material::{
        lineage::{self as run_lineage, EventBinding},
        run::{self as run_material, MaterialOrigin},
    },
    sequence::analysis::{self, Primer},
    verification,
};
use base64::{Engine as _, engine::general_purpose::STANDARD};
use collagene_model::{Alphabet, Error, Molecule};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};
pub const ENGINE_VERSION: &str = "material-derivation-5";
pub const WORKSPACE_ATTACHMENT_KIND: &str = "material-derivation-v5";
pub const MEDIA_TYPE: &str = "application/vnd.collagene.material-derivation";
pub const EVENT_OPERATION: &str = "Target-free native material derivation v5";
pub const MAX_RECORD_BYTES: usize = 8 * 1024 * 1024;
pub const MAX_CAPSULE_BYTES: usize = 32 * 1024 * 1024;
pub const MAX_RECORDS: usize = 128;
pub const MAX_SOURCES: usize = 2;
const MAGIC: &[u8; 8] = b"MDRV\0\0\0\x05";
const HEADER: usize = 44;
const CONTEXT_KEY: &str = "source:material_derivation:context";
macro_rules! dto{($($t:tt)*)=>{#[derive(Debug,Clone,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
mod ancestry;
pub mod annealing;
pub mod digest;
pub mod end_repair;
pub mod overlap_extension;
pub mod pcr;
pub mod terminal;

dto! {#[derive(Copy,PartialEq,Eq,PartialOrd,Ord)]#[serde(rename_all="snake_case")]pub enum SourceRole{Template,TopOligo,BottomOligo}}
dto! {pub struct BoundSource{pub role:SourceRole,pub event:EventBinding,pub origin:SourceOrigin}}
dto! {#[derive(PartialEq,Eq)]pub struct SourceBinding{pub role:SourceRole,pub event:EventBinding}}
dto! {pub struct SourceSnapshot{pub role:SourceRole,pub molecule:Molecule}}
dto! {pub enum SourceOrigin{OriginalReference,Saved(MaterialOrigin)}}
dto! {pub struct PcrSelection{pub selected_index:usize,pub expected_preview_revision:String,#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]pub geometry:analysis::PcrProduct}}
dto! {pub struct MismatchPcrSelection{pub selected_index:usize,pub expected_preview_revision:String,#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]pub geometry:crate::cloning::pcr::mutagenesis::MutantPcrProduct}}
dto! {pub struct OverlapExtensionSelection{pub expected_preview_revision:String,pub expected_product_revision:String}}
dto! {#[derive(Copy,PartialEq,Eq)]pub enum TerminalAction{TreatEnds{treatment:chemistry::EndTreatment},Tail3Prime{base:crate::cloning::ligation::tailing::TailingBase}}}
dto! {pub struct AnnealSelection{pub expected_preview_revision:String,pub candidate_id:String,pub expected_product_revision:String}}
dto! {
    pub enum Derivation {
        ExactPcr {
            forward: Primer,
            reverse: Primer,
            template_strand: chemistry::Strand,
            selection: PcrSelection,
        },
        CatalogDigest {
            selection: CatalogSelection,
            options: InspectionOptions,
            #[serde(deserialize_with = "crate::encoding::required_option")]
            #[cfg_attr(feature = "json-schema", schemars(required))]
            conditions: Option<crate::cloning::restriction::conditions::ReactionConditions>,
            fragment_index: usize,
        },
        MismatchPcr {
            forward: Primer,
            reverse: Primer,
            options: crate::cloning::pcr::mutagenesis::PcrOptions,
            selection: MismatchPcrSelection,
        },
        EndRepair { method: chemistry::EndRepairMethod },
        OverlapExtensionPcr {
            preparation: crate::cloning::overlap::extension::pcr::Request,
            selection: OverlapExtensionSelection,
        },
        TerminalProcessing { action: TerminalAction },
        AnnealedOligos {
            options: chemistry::AnnealOptions,
            selection: AnnealSelection,
        },
    }
}
dto! {pub struct Request{pub sources:Vec<BoundSource>,pub derivation:Derivation,pub name:String,pub save_primers:bool}}
dto! {pub struct Record{pub schema_version:u32,pub engine_version:String,pub scientific_engine:String,pub source_workspace_fingerprint:String,pub context_revision:String,pub request:Request,pub product:EventBinding,pub product_state_revision:String,pub incomplete:bool}}
mod record_bytes {
    use super::*;
    use serde::{
        Deserializer, Serializer,
        de::{self, Visitor},
    };
    use std::fmt;
    pub fn serialize<S: Serializer>(bytes: &[u8], s: S) -> Result<S::Ok, S::Error> {
        if bytes.len() > MAX_RECORD_BYTES {
            return Err(serde::ser::Error::custom("derived record byte bound"));
        }
        if s.is_human_readable() {
            s.serialize_str(&STANDARD.encode(bytes))
        } else {
            s.serialize_bytes(bytes)
        }
    }
    pub fn deserialize<'de, D: Deserializer<'de>>(d: D) -> Result<Vec<u8>, D::Error> {
        struct Bytes {
            human: bool,
        }
        impl<'de> Visitor<'de> for Bytes {
            type Value = Vec<u8>;
            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                f.write_str("bounded base64 or binary native derivation record")
            }
            fn visit_str<E: de::Error>(self, s: &str) -> Result<Vec<u8>, E> {
                if !self.human || s.len() > MAX_RECORD_BYTES.div_ceil(3) * 4 {
                    return Err(E::custom("record byte bound"));
                }
                let bytes = STANDARD.decode(s).map_err(E::custom)?;
                if bytes.len() > MAX_RECORD_BYTES {
                    return Err(E::custom("record byte bound"));
                }
                Ok(bytes)
            }
            fn visit_bytes<E: de::Error>(self, bytes: &[u8]) -> Result<Vec<u8>, E> {
                if bytes.len() > MAX_RECORD_BYTES {
                    return Err(E::custom("record byte bound"));
                }
                Ok(bytes.to_vec())
            }
            fn visit_byte_buf<E: de::Error>(self, bytes: Vec<u8>) -> Result<Vec<u8>, E> {
                if bytes.len() > MAX_RECORD_BYTES {
                    return Err(E::custom("record byte bound"));
                }
                Ok(bytes)
            }
        }
        let human = d.is_human_readable();
        if human {
            d.deserialize_str(Bytes { human })
        } else {
            d.deserialize_bytes(Bytes { human })
        }
    }
}
dto! {pub struct Capsule{pub schema_version:u16,pub product:EventBinding,#[serde(with="record_bytes")]#[cfg_attr(feature="json-schema",schemars(with="String"))]pub record_bytes:Vec<u8>,pub source_snapshots:Vec<SourceSnapshot>}}
dto! {
    pub struct Receipt {
        pub engine_version: String,
        pub context_revision: String,
        pub record_hash: String,
        pub primer_attachment: Option<String>,
        pub sources: Vec<SourceBinding>,
        pub product: EventBinding,
        pub incomplete: bool,
        pub findings: Vec<String>,
        pub findings_total: usize,
        pub findings_truncated: bool,
        pub duplex: chemistry::DuplexSummary,
        pub modifications: ModificationSummary,
        /// Bounded originating-report provenance from this fresh execution;
        /// serialized receipts remain documentary, never reusable authority.
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, crate::cloning::restriction::conditions::source::SourceReport, 32>")]
        pub reaction_conditions: Vec<crate::cloning::restriction::conditions::source::SourceReport>,
    }
}
pub struct PreparedInstallation {
    pub workspace: Workspace,
    pub receipt: Receipt,
}
dto! {pub struct Installed{pub archive:Vec<u8>,pub receipt:Receipt}}
/// Not serializable, not publicly constructible, and never reconstructed from summaries.
#[derive(Debug)]
pub struct FreshMaterial {
    reference: Molecule,
    tracked: TrackedDuplex,
    incomplete: bool,
    findings: Vec<String>,
    findings_total: usize,
    findings_truncated: bool,
    reaction_conditions: Vec<crate::cloning::restriction::conditions::source::SourceReport>,
}
impl FreshMaterial {
    fn from_run_source(fresh: run_material::FreshMaterial) -> Result<Self, Error> {
        let reaction_conditions = reaction_conditions::collect(fresh.reaction_conditions())?;
        let (reference, tracked, incomplete, findings, findings_total, findings_truncated) =
            fresh.into_parts();
        Ok(Self {
            reference,
            tracked,
            incomplete,
            findings,
            findings_total,
            findings_truncated,
            reaction_conditions,
        })
    }

    pub fn reaction_conditions(
        &self,
    ) -> &[crate::cloning::restriction::conditions::source::SourceReport] {
        &self.reaction_conditions
    }
    pub fn reference(&self) -> &Molecule {
        &self.reference
    }
    pub fn tracked(&self) -> &TrackedDuplex {
        &self.tracked
    }
    pub fn incomplete(&self) -> bool {
        self.incomplete
    }
    pub fn findings(&self) -> &[String] {
        &self.findings
    }
    pub fn findings_total(&self) -> usize {
        self.findings_total
    }
    pub fn findings_truncated(&self) -> bool {
        self.findings_truncated
    }
    pub fn into_parts(self) -> (Molecule, TrackedDuplex, bool, Vec<String>, usize, bool) {
        (
            self.reference,
            self.tracked,
            self.incomplete,
            self.findings,
            self.findings_total,
            self.findings_truncated,
        )
    }
}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("native material derivation: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
struct Bound {
    bytes: Vec<u8>,
    limit: usize,
}
impl Write for Bound {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if self
            .bytes
            .len()
            .checked_add(b.len())
            .is_none_or(|n| n > self.limit)
        {
            return Err(std::io::Error::other("derivation record bound"));
        }
        self.bytes.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn encoded(value: &impl Serialize, limit: usize) -> Result<Vec<u8>, Error> {
    let mut out = Bound {
        bytes: Vec::new(),
        limit,
    };
    serde_json::to_writer(&mut out, value).map_err(|_| Error::TooLarge)?;
    Ok(out.bytes)
}
fn name(s: &str) -> Result<(), Error> {
    if s.trim().is_empty() || s.trim() != s || s.len() > 256 || s.chars().any(char::is_control) {
        return Err(bad("declare a bounded nonempty product name"));
    }
    Ok(())
}
impl SourceOrigin {
    pub fn validate_shallow_for(&self, event: &EventBinding) -> Result<(), Error> {
        run_material::event_shape(event)?;
        if let Self::Saved(origin) = self {
            if origin.event() != event {
                return Err(bad("origin differs from exact selected source event"));
            }
        }
        Ok(())
    }
}
impl Request {
    pub fn template(&self) -> Result<&BoundSource, Error> {
        if self.sources.len() != 1 || self.sources[0].role != SourceRole::Template {
            return Err(bad(
                "this derivation requires exactly the Template source role",
            ));
        }
        Ok(&self.sources[0])
    }
    pub fn source_bindings(&self) -> Vec<SourceBinding> {
        self.sources
            .iter()
            .map(|s| SourceBinding {
                role: s.role,
                event: s.event.clone(),
            })
            .collect()
    }
    fn validate_sources(&self) -> Result<(), Error> {
        let expected: &[SourceRole] =
            if matches!(self.derivation, Derivation::AnnealedOligos { .. }) {
                &[SourceRole::TopOligo, SourceRole::BottomOligo]
            } else {
                &[SourceRole::Template]
            };
        if self.sources.len() != expected.len() {
            return Err(bad(
                "native source roles/count differ from the declared derivation",
            ));
        }
        for (source, role) in self.sources.iter().zip(expected) {
            if source.role != *role {
                return Err(bad(
                    "native source roles must be unique and in canonical method order",
                ));
            }
            source.origin.validate_shallow_for(&source.event)?;
            if matches!(self.derivation, Derivation::AnnealedOligos { .. })
                && !matches!(source.origin, SourceOrigin::OriginalReference)
            {
                return Err(bad(
                    "Annealing requires independently declared standard oligos; saved duplexes are not denatured or reinterpreted as single strands",
                ));
            }
        }
        Ok(())
    }

    /// One level only. The common mixed walker validates saved descendants under
    /// one ledger before calling this method; never recurse through serde here.
    pub(crate) fn validate_shallow(&self) -> Result<(), Error> {
        name(&self.name)?;
        self.validate_sources()?;
        match &self.derivation {
            Derivation::AnnealedOligos { options, selection } => {
                if self.save_primers {
                    return Err(bad(
                        "Annealing uses the declared oligo inputs and creates no new primer order",
                    ));
                }
                annealing::anneal_options_shape(options)?;
                if !hex(&selection.expected_preview_revision)
                    || !hex(&selection.expected_product_revision)
                    || !annealing::candidate_id_shape(&selection.candidate_id)
                {
                    return Err(bad("invalid bound annealing candidate selection"));
                }
            }

            Derivation::ExactPcr {
                forward,
                reverse,
                template_strand,
                selection,
            } => {
                for p in [forward, reverse] {
                    if p.annealing.len() > 1000 || p.tail.len() > 1000 {
                        return Err(Error::TooLarge);
                    }
                    crate::cloning::pcr::outputs::primer(p)?;
                    if p.ordered_sequence().bytes().any(|b| !b"ACGT".contains(&b)) {
                        return Err(bad(
                            "reusable exact PCR requires concrete newly synthesized primer bases",
                        ));
                    }
                }
                if selection.geometry.template_strand != *template_strand
                    || selection.selected_index >= 10000
                    || !hex(&selection.expected_preview_revision)
                    || selection.geometry.sequence.is_empty()
                    || selection.geometry.sequence.len() > collagene_model::MAX_SEQUENCE_LEN
                    || selection.geometry.sequence.alphabet() != Alphabet::Dna
                    || selection
                        .geometry
                        .sequence
                        .as_str()
                        .bytes()
                        .any(|b| !b"ACGT".contains(&b))
                    || selection.geometry.forward_primer.len() > 256
                    || selection.geometry.reverse_primer.len() > 256
                {
                    return Err(bad("invalid source-bound PCR selection"));
                }
            }
            Derivation::MismatchPcr {
                forward,
                reverse,
                options,
                selection,
            } => {
                crate::cloning::assembly::recipe::mismatch_pcr_shape(forward, reverse, options)?;
                if selection.selected_index >= 10000
                    || !hex(&selection.expected_preview_revision)
                    || selection.geometry.sequence.is_empty()
                    || selection.geometry.sequence.len() > collagene_model::MAX_SEQUENCE_LEN
                    || selection.geometry.sequence.alphabet() != Alphabet::Dna
                    || selection
                        .geometry
                        .sequence
                        .as_str()
                        .bytes()
                        .any(|b| !b"ACGT".contains(&b))
                    || selection.geometry.substitutions.len() > 400
                    || selection.geometry.forward_binding.mismatches.len() > 200
                    || selection.geometry.reverse_binding.mismatches.len() > 200
                {
                    return Err(bad(
                        "invalid source-bound intentional-mismatch PCR selection",
                    ));
                }
            }
            Derivation::OverlapExtensionPcr {
                preparation,
                selection,
            } => {
                preparation.validate_shape()?;
                if !hex(&selection.expected_product_revision)
                    || selection.expected_preview_revision
                        != overlap_extension::overlap_extension_context(
                            &self.template()?.event,
                            &self.template()?.origin,
                            preparation,
                        )?
                {
                    return Err(bad(
                        "stale or invalid exact-source overlap-extension selection",
                    ));
                }
                if matches!(self.template()?.origin, SourceOrigin::OriginalReference)
                    && [&preparation.left, &preparation.right]
                        .iter()
                        .any(|pair| pair.template_strand != chemistry::Strand::Top)
                {
                    return Err(bad(
                        "Original-reference OE PCR requires Top source strands; Bottom requires native material replay",
                    ));
                }
            }
            Derivation::TerminalProcessing { .. } => {
                if self.save_primers {
                    return Err(bad(
                        "Terminal processing declares no new oligos; historical source primers remain ancestry",
                    ));
                }
            }
            Derivation::EndRepair { .. } => {
                if self.save_primers {
                    return Err(bad(
                        "End repair creates no new primer order; historical primers remain source ancestry.",
                    ));
                }
            }
            Derivation::CatalogDigest {
                selection,
                options,
                conditions,
                fragment_index,
            } => {
                if self.save_primers {
                    return Err(bad(
                        "Catalog digest cannot create a new PCR primer order; existing upstream primers remain source history.",
                    ));
                }
                selection.validate()?;
                options.validate()?;
                if let Some(conditions) = conditions {
                    conditions.validate_for(&catalog_digest::builtin_catalog()?, selection)?;
                }
                if *fragment_index >= options.max_fragments {
                    return Err(bad("fragment index exceeds declared bounds"));
                }
            }
        }
        Ok(())
    }
}
fn context(fingerprint: &str, request: &Request, science: &str) -> Result<String, Error> {
    crate::encoding::bounded_json::hash(&(ENGINE_VERSION, science, fingerprint, request))
}
impl Record {
    pub fn validate_documentary(&self) -> Result<(), Error> {
        crate::material::capsule::derived_record_preflight(self)
    }
    pub(crate) fn validate_shallow(&self) -> Result<(), Error> {
        if self.schema_version != 5
            || self.engine_version != ENGINE_VERSION
            || self.scientific_engine.is_empty()
            || self.scientific_engine.len() > 128
            || !self.scientific_engine.is_ascii()
            || self.scientific_engine.bytes().any(|b| b.is_ascii_control())
            || !hex(&self.source_workspace_fingerprint)
            || !hex(&self.context_revision)
            || !hex(&self.product_state_revision)
        {
            return Err(bad("unsupported record schema/engine or fingerprint"));
        }
        self.request.validate_shallow()?;
        run_material::event_shape(&self.product)?;
        if self
            .request
            .sources
            .iter()
            .any(|s| self.product.event_id == s.event.event_id)
            || context(
                &self.source_workspace_fingerprint,
                &self.request,
                &self.scientific_engine,
            )? != self.context_revision
        {
            return Err(bad("record source/product/context differs"));
        }
        Ok(())
    }
    pub fn source_bindings(&self) -> Result<Vec<EventBinding>, Error> {
        self.validate_documentary()?;
        let mut values = BTreeMap::new();
        let mut add = |e: &EventBinding| {
            values.insert(
                (
                    e.document_id.clone(),
                    e.event_id.clone(),
                    e.revision.clone(),
                ),
                e.clone(),
            );
        };
        for source in &self.request.sources {
            add(&source.event);
        }
        // Collect callback results separately so only one closure mutably borrows the map.
        let mut runs = Vec::new();
        let mut derived = Vec::new();
        let mut graphs = Vec::new();
        crate::material::capsule::derived_ancestry_walk(
            self,
            |origin, _| {
                runs.push(origin.target_event.clone());
                runs.push(origin.selection.event.clone());
                runs.extend(origin.inputs.iter().map(|b| b.event.clone()));
                runs.extend(origin.products.iter().map(|p| EventBinding {
                    document_id: p.document_id.clone(),
                    event_id: p.event_id.clone(),
                    revision: p.revision.clone(),
                }));
                Ok(())
            },
            |r, _| {
                derived.extend(r.request.sources.iter().map(|s| s.event.clone()));
                if r.product != self.product {
                    derived.push(r.product.clone())
                }
                Ok(())
            },
            |graph, _| {
                graphs.extend(graph.request.bindings.iter().map(|b| b.event.clone()));
                Ok(())
            },
        )?;
        for e in runs.iter().chain(&derived).chain(&graphs) {
            add(e)
        }
        Ok(values.into_values().collect())
    }
    pub fn required_run_ids(&self) -> Result<Vec<String>, Error> {
        let mut ids = BTreeSet::new();
        crate::material::capsule::derived_ancestry_walk(
            self,
            |origin, _| {
                ids.insert(origin.selection.run_id.clone());
                Ok(())
            },
            |_, _| Ok(()),
            |_, _| Ok(()),
        )?;
        Ok(ids.into_iter().collect())
    }
    pub fn required_derivation_products(&self) -> Result<Vec<EventBinding>, Error> {
        let mut ids = BTreeMap::new();
        crate::material::capsule::derived_ancestry_walk(
            self,
            |_, _| Ok(()),
            |r, _| {
                if r.product != self.product {
                    ids.insert(r.product.event_id.clone(), r.product.clone());
                }
                Ok(())
            },
            |_, _| Ok(()),
        )?;
        Ok(ids.into_values().collect())
    }
    /// Exact saved TOPO graph dependencies retained by this material ancestry.
    pub fn required_graph_records(
        &self,
    ) -> Result<Vec<crate::cloning::topo::graph_history::Selection>, Error> {
        let mut records = BTreeMap::new();
        crate::material::capsule::derived_ancestry_walk(
            self,
            |_, _| Ok(()),
            |_, _| Ok(()),
            |graph, _| {
                let bytes = crate::cloning::topo::graph_history::save_record(graph)?;
                let hash = blake3::hash(&bytes).to_hex().to_string();
                records.insert(
                    hash.clone(),
                    crate::cloning::topo::graph_history::Selection {
                        record_hash: hash,
                        state_revision: graph.state_revision.clone(),
                    },
                );
                Ok(())
            },
        )?;
        Ok(records.into_values().collect())
    }
}
fn save_record_shallow(record: &Record) -> Result<Vec<u8>, Error> {
    record.validate_shallow()?;
    let payload = encoded(record, MAX_RECORD_BYTES - HEADER)?;
    let mut bytes = Vec::with_capacity(HEADER + payload.len());
    bytes.extend_from_slice(MAGIC);
    bytes.extend_from_slice(&(payload.len() as u32).to_le_bytes());
    bytes.extend_from_slice(blake3::hash(&payload).as_bytes());
    bytes.extend_from_slice(&payload);
    Ok(bytes)
}
pub fn save_record(record: &Record) -> Result<Vec<u8>, Error> {
    record.validate_documentary()?;
    save_record_shallow(record)
}
/// Header/size/checksum admission before decoding exactly ONE JSON level. Saved
/// Derived origins contain opaque bytes, never recursively nested Record DTOs.
pub(crate) fn decode_record_resources(bytes: &[u8]) -> Result<Record, Error> {
    if bytes.len() < HEADER || bytes.len() > MAX_RECORD_BYTES || &bytes[..8] != MAGIC {
        return Err(bad("invalid bounded native record header"));
    }
    let length = u32::from_le_bytes(bytes[8..12].try_into().unwrap()) as usize;
    if length != bytes.len() - HEADER || blake3::hash(&bytes[HEADER..]).as_bytes() != &bytes[12..44]
    {
        return Err(bad("native record size/checksum differs"));
    }
    let record: Record = serde_json::from_slice(&bytes[HEADER..])
        .map_err(|_| bad("invalid native derivation JSON"))?;
    if record.schema_version != 5 || record.request.sources.len() > MAX_SOURCES {
        return Err(bad("unsupported or oversized native derivation record"));
    }
    Ok(record)
}
/// Current-frame admission followed by exact documentary semantics.
pub(crate) fn load_record_shallow(bytes: &[u8]) -> Result<Record, Error> {
    let record = decode_record_resources(bytes)?;
    record.validate_shallow()?;
    if save_record_shallow(&record)? != bytes {
        return Err(bad("record is not canonical"));
    }
    Ok(record)
}
pub fn load_record(bytes: &[u8]) -> Result<Record, Error> {
    let record = load_record_shallow(bytes)?;
    record.validate_documentary()?;
    Ok(record)
}
fn original_reference_admission(source: &Molecule) -> Result<(), Error> {
    collagene_model::bounds::model_preflight(&[source])?;
    source.validate()?;
    crate::formats::fidelity::reference_export_guard(source)?;
    if source.sequence.alphabet() != Alphabet::Dna {
        return Err(bad("PCR/digest material requires DNA"));
    }
    Ok(())
}
fn original_source(source: &Molecule) -> Result<TrackedDuplex, Error> {
    original_reference_admission(source)?;
    // Plain originals are an explicit intact-reference hypothesis, never an
    // imported chemical state. Unknown source modifications are not set to zero.
    TrackedDuplex::from_unknown_duplex(&chemistry::Duplex::from_molecule(source)?)
}
fn limitations(
    findings: Vec<String>,
    total: usize,
    already_truncated: bool,
) -> (Vec<String>, usize, bool) {
    let count = total.max(findings.len());
    let mut truncated = already_truncated || findings.len() > 128;
    let values = findings
        .into_iter()
        .take(128)
        .map(|s| {
            if s.chars().count() > 2048 {
                truncated = true;
                s.chars().take(2048).collect()
            } else {
                s
            }
        })
        .collect();
    (values, count, truncated)
}
/// Only native, freshly executed sources may reach this adapter. Root's generic
/// material resolver preserves both saved Run and Derived ancestry.
pub(crate) fn fresh_source(
    source: &Molecule,
    origin: &SourceOrigin,
) -> Result<FreshMaterial, Error> {
    collagene_model::bounds::model_preflight(&[source])?;
    source.validate()?;
    let (tracked,incomplete,findings,total,truncated,reaction_conditions)=match origin {
        SourceOrigin::OriginalReference=>(original_source(source)?,false,vec!["Original source is an explicit intact DNA reference hypothesis with unknown source modifications and terminal groups; no measured stock chemistry or amplification success is inferred.".into()],1,false,vec![]),
        SourceOrigin::Saved(origin)=>{
            let fresh=origin.replay(source)?;
            let evidence=reaction_conditions::collect(fresh.reaction_conditions())?;
            let (_,tracked,incomplete,findings,total,truncated)=fresh.into_parts();
            (tracked,incomplete,findings,total,truncated,evidence)
        },
    };
    let (findings, findings_total, findings_truncated) = limitations(findings, total, truncated);
    Ok(FreshMaterial {
        reaction_conditions,
        reference: source.clone(),
        tracked,
        incomplete,
        findings,
        findings_total,
        findings_truncated,
    })
}
fn request_source(source: &Molecule, request: &Request) -> Result<(), Error> {
    let template = request.template()?;
    collagene_model::bounds::model_preflight(&[source])?;
    source.validate()?;
    request.validate_shallow()?;
    if source.revision() != template.event.revision {
        return Err(bad("source snapshot differs from exact request event"));
    }
    if let Derivation::MismatchPcr {
        forward,
        reverse,
        options,
        ..
    } = &request.derivation
    {
        if crate::cloning::assembly::recipe::mismatch_pcr_work(
            source.sequence.len(),
            forward,
            reverse,
            options,
        )? > 200_000_000
        {
            return Err(Error::TooLarge);
        }
    }
    if let Derivation::OverlapExtensionPcr { preparation, .. } = &request.derivation {
        crate::cloning::overlap::extension::execution::preflight(source, preparation)?;
    }
    if let SourceOrigin::OriginalReference = &template.origin {
        original_source(source)?;
    }
    Ok(())
}
pub fn preflight(w: &Workspace, request: &Request) -> Result<(), Error> {
    request.validate_shallow()?;
    w.validate()?;
    crate::material::capsule::validate_derived_request_origins(w, request)?;
    for source in &request.sources {
        let molecule = run_material::resolve(w, &source.event)?;
        if matches!(request.derivation, Derivation::AnnealedOligos { .. }) {
            declared_oligo_source(molecule)?;
        } else {
            request_source(molecule, request)?;
        }
        let doc = w
            .documents
            .iter()
            .find(|d| d.id == source.event.document_id)
            .ok_or_else(|| bad("source document absent"))?;
        if doc.head != source.event.revision
            || doc.head_event.as_deref() != Some(source.event.event_id.as_str())
        {
            return Err(bad("source head changed; recapture the selected document"));
        }
    }
    if w.attachments.len() >= crate::documents::workspace::MAX_ATTACHMENTS
        || w.attachments
            .values()
            .filter(|a| a.kind == WORKSPACE_ATTACHMENT_KIND)
            .count()
            >= MAX_RECORDS
    {
        return Err(Error::TooLarge);
    }
    Ok(())
}
pub fn revision(w: &Workspace, request: &Request) -> Result<String, Error> {
    preflight(w, request)?;
    context(
        &run_lineage::workspace_fingerprint(w)?,
        request,
        verification::ENGINE_VERSION,
    )
}
fn fresh_template(
    source: &Molecule,
    request: &Request,
    revision: &str,
    capture: &mut Option<ancestry::Capture>,
    replay: &mut verification::replay_budget::ReplayBudget<'_>,
) -> Result<FreshMaterial, Error> {
    let template = request.template()?;
    request_source(source, request)?;
    let mut inherited = match capture {
        Some(capture) => capture.source(template.role, source, &template.origin, replay)?,
        None if replay.is_qc_job() => match &template.origin {
            SourceOrigin::OriginalReference => fresh_source(source, &template.origin)?,
            SourceOrigin::Saved(origin) => {
                let (fresh, _) = origin.replay_scoped(source, replay, None)?;
                FreshMaterial::from_run_source(fresh)?
            }
        },
        None => fresh_source(source, &template.origin)?,
    };
    let (mut reference, tracked, local) = match &request.derivation {
        Derivation::AnnealedOligos { .. } => {
            return Err(bad("annealing requires the plural oligo source path"));
        }
        Derivation::ExactPcr {
            forward,
            reverse,
            template_strand,
            selection,
        } => {
            let previews = if replay.is_qc_job() {
                match &template.origin {
                    SourceOrigin::OriginalReference => {
                        verification::pcr::simulate_pcr_documents_with_budget(
                            source,
                            forward,
                            reverse,
                            *template_strand,
                            &mut replay.pcr,
                        )?
                    }
                    SourceOrigin::Saved(origin) => {
                        crate::cloning::pcr::material::simulate_from_fresh(
                            source,
                            origin,
                            forward,
                            reverse,
                            *template_strand,
                            &inherited.tracked,
                            inherited.incomplete,
                            &inherited.findings,
                            inherited.findings_total,
                            inherited.findings_truncated,
                            Some(&mut replay.pcr),
                        )?
                    }
                }
            } else {
                let previews = match &template.origin {
                    SourceOrigin::OriginalReference => verification::pcr::simulate_pcr_documents(
                        source,
                        forward,
                        reverse,
                        *template_strand,
                    )?,
                    SourceOrigin::Saved(origin) => {
                        crate::cloning::pcr::material::simulate_from_fresh(
                            source,
                            origin,
                            forward,
                            reverse,
                            *template_strand,
                            &inherited.tracked,
                            inherited.incomplete,
                            &inherited.findings,
                            inherited.findings_total,
                            inherited.findings_truncated,
                            None,
                        )?
                    }
                };
                previews
            };
            if previews.len() != 1 || selection.selected_index != 0 {
                return Err(bad(
                    "initial reusable PCR profile requires exactly one complete geometry; multiple products remain reference previews, never inferred isolation",
                ));
            }
            let preview = &previews[0];
            if preview.molecule.revision() != selection.expected_preview_revision
                || crate::encoding::bounded_json::hash(&preview.geometry)?
                    != crate::encoding::bounded_json::hash(&selection.geometry)?
            {
                return Err(bad("selected PCR preview/geometry changed"));
            }
            let frame = crate::cloning::pcr::source::Template::fresh(
                source,
                &inherited.tracked,
                *template_strand,
            )?;
            let geometries = if replay.is_qc_job() {
                frame
                    .geometries_with_budget(forward, reverse, &mut replay.pcr)?
                    .0
            } else {
                frame.geometries(forward, reverse)?
            };
            if geometries.len() != 1 {
                return Err(bad("actual exact PCR geometry is not unique"));
            }
            let geometry = &geometries[0];
            let tracked = TrackedDuplex::exact_pcr(&geometry.product.sequence)?;
            let mut reference = tracked.duplex.molecule(&request.name)?;
            let mut findings = Vec::new();
            let receipt = verification::pcr::project_pcr_frame_annotations(
                frame,
                geometry,
                &mut reference,
                &mut findings,
            )?;
            if let Some(capture) = capture {
                capture.record(
                    0,
                    receipt,
                    verification::ancestry::ProjectionKind::TemplateSynthesis,
                    vec![],
                    verification::ancestry::IntroductionKind::PrimerTail,
                );
            }
            findings.push("Unique exact PCR produces newly synthesized unmodified blunt duplex strands with 5′OH/3′OH in the existing ideal-PCR model. Primer chemistry, yields, cleanup, amplification efficiency and experimental isolation are not measured.".into());
            (reference, tracked, findings)
        }
        Derivation::MismatchPcr {
            forward,
            reverse,
            options,
            selection,
        } => {
            let preview = pcr::mismatch_from_fresh(
                source,
                &template.origin,
                &inherited,
                forward,
                reverse,
                options,
            )?;
            if preview.preview.report.products.len() != 1 || selection.selected_index != 0 {
                return Err(bad(
                    "Reusable intentional-mismatch PCR requires exactly one reported native geometry; competing products remain previews, not inferred isolation.",
                ));
            }
            let geometry = &preview.preview.report.products[0];
            if !geometry.certain || geometry.overlapping_domains {
                return Err(bad(
                    "Ambiguous/overlapping mismatch PCR geometry remains diagnostic, not established amplified material.",
                ));
            }
            if preview.preview.documents[0].revision() != selection.expected_preview_revision
                || crate::encoding::bounded_json::hash(geometry)?
                    != crate::encoding::bounded_json::hash(&selection.geometry)?
            {
                return Err(bad("selected mismatch-PCR preview/geometry changed"));
            }
            inherited.incomplete |= !preview.preview.report.enumeration_complete;
            let tracked = TrackedDuplex::exact_pcr(&geometry.sequence)?;
            let mut reference = tracked.duplex.molecule(&request.name)?;
            let mut findings = preview.preview.report.warnings.clone();
            let receipt = verification::pcr::project_mutant_pcr(
                source,
                geometry,
                &mut reference,
                &mut findings,
            )?;
            if let Some(capture) = capture {
                let changes = geometry
                    .substitutions
                    .iter()
                    .map(|change| verification::ancestry::MappedChange {
                        source_position: change.source_position,
                        destination_position: change.product_position,
                        before: change.template_base,
                        after: change.product_base,
                    })
                    .collect();
                capture.record(
                    0,
                    receipt,
                    verification::ancestry::ProjectionKind::EditedTemplate,
                    changes,
                    verification::ancestry::IntroductionKind::PrimerTail,
                );
            }
            findings.push(format!("Intentional-mismatch PCR retains {} primer-encoded substitutions, both complete ordered oligos and the declared 3′/mismatch/search policy. The selected native geometry supplies newly synthesized unmodified blunt strands with 5′OH/3′OH; amplification efficiency, isolation and target success are not established.",geometry.substitutions.len()));
            if !preview.preview.report.enumeration_complete {
                findings.push("Alternative mismatch-PCR search is incomplete. The selected singleton path remains conditional material and its limitations are not upgraded.".into());
            }
            (reference, tracked, findings)
        }
        Derivation::OverlapExtensionPcr {
            preparation,
            selection,
        } => {
            let executed = crate::cloning::overlap::extension::execution::execute_scoped(
                source,
                Some(&inherited.tracked),
                preparation,
                capture.is_some(),
                replay,
            )?;
            let mut product=executed.product.ok_or_else(||bad("Reusable OE material requires one actual geometry in each source PCR and outer PCR; no product is inferred from a desired sequence"))?;
            if product.revision() != selection.expected_product_revision {
                return Err(bad(
                    "Selected outer-PCR product changed; inspect the source-backed OE preparation again",
                ));
            }
            product.name = request.name.clone();
            if let Some(capture) = capture {
                capture.record_oe(executed.native);
            }
            inherited.incomplete |= executed.incomplete;
            let tracked = TrackedDuplex::exact_pcr(&product.sequence)?;
            let mut findings = executed.findings;
            findings.push("Three literal primer pairs produce the actual outer-PCR product. Newly synthesized unmodified blunt strands carry 5′OH/3′OH under the ideal-PCR model. Source/terminal uncertainty remains inherited; amplification efficiency and experimental isolation are not established, and no intended target was checked.".into());
            (product, tracked, findings)
        }
        Derivation::TerminalProcessing { action } => {
            let result = terminal::terminal_from_fresh(source, &inherited, action, &request.name)?;
            if let Some(capture) = capture {
                capture.record_top_copy(source, &result.0)?;
            }
            result
        }
        Derivation::EndRepair { method } => {
            let repair = chemistry::repair_ends(&inherited.tracked.duplex, *method)?;
            let tracked = inherited.tracked.repair(*method)?;
            if tracked.duplex != repair.duplex {
                return Err(bad("native repair geometry and tracked mappings differ"));
            }
            let mut findings = repair.warnings.clone();
            let projected = verification::terminal::repaired_document(
                source,
                &repair,
                request.name.clone(),
                &mut findings,
            )?;
            let reference = projected.molecule;
            if let Some(capture) = capture {
                for receipt in projected.projections {
                    capture.record(
                        0,
                        receipt,
                        verification::ancestry::ProjectionKind::ReferenceCorrespondence,
                        vec![],
                        verification::ancestry::IntroductionKind::NewlySynthesizedUnmappedCalls,
                    );
                }
            }
            findings.push("Direct end repair uses actual source strands and template-supported fill/removal only. Surviving positional modifications, original 5′ groups and internal nicks remain; no phosphorylation, nick sealing, supplier-vector chemistry or target success is inferred.".into());
            (reference, tracked, findings)
        }
        Derivation::CatalogDigest {
            selection,
            options,
            conditions,
            fragment_index,
        } => {
            let catalog = catalog_digest::builtin_catalog()?;
            selection.validate_for(&catalog)?;
            if let Some(options) = conditions {
                let evidence =
                    crate::cloning::restriction::conditions::source::SourceReport::assess(
                        &catalog,
                        source,
                        &inherited.tracked,
                        selection,
                        options,
                    )?;
                inherited.incomplete |= !evidence.satisfies_request();
                inherited.reaction_conditions = reaction_conditions::collect(
                    inherited
                        .reaction_conditions
                        .iter()
                        .chain(std::iter::once(&evidence)),
                )?;
            }
            let enzymes = selection.resolved_enzymes(&catalog)?;
            let role = crate::cloning::restriction::history::source_role(&template.event)?;
            let mut selected = None;
            let inspection = restriction_inspection::inspect_with_digest_for_source(
                &role,
                source,
                &enzymes,
                options,
                &inherited.tracked.duplex.summary()?,
                || {
                    let mut result = Vec::new();
                    let fragments = if replay.is_qc_job() {
                        modified_duplex::digest_with_budget(
                            &inherited.tracked,
                            &enzymes,
                            replay.cleavage,
                        )?
                    } else {
                        modified_duplex::digest(&inherited.tracked, &enzymes)?
                    };
                    for (index, fragment) in fragments.into_iter().enumerate() {
                        result.push(chemistry::DigestFragment {
                            duplex: fragment.duplex.duplex.clone(),
                            top_start: fragment.top_start,
                            bottom_start: fragment.bottom_start,
                        });
                        if index == *fragment_index {
                            selected = Some(fragment)
                        }
                    }
                    Ok(result)
                },
            )?;
            if !inspection.products_determined {
                return Err(bad(
                    "bounded actual cleavage does not determine a reusable fragment",
                ));
            }
            let fragment = selected.ok_or_else(|| bad("selected actual digest fragment absent"))?;
            let mut reference = fragment.duplex.duplex.molecule(&request.name)?;
            let geometry = chemistry::DigestFragment {
                duplex: fragment.duplex.duplex.clone(),
                top_start: fragment.top_start,
                bottom_start: fragment.bottom_start,
            };
            let mut findings = inspection.diagnostics;
            if conditions.is_some() {
                if let Some(evidence) = inherited.reaction_conditions.iter().find(|e| {
                    e.source_revision == source.revision()
                        && e.report.selection == *selection
                        && conditions.as_ref() == Some(&e.report.options)
                }) {
                    findings.push(format!("Retained source-bound documentary reaction constraint: {:?}. Star activity and incubation-duration compatibility remain unassessed; sequence product geometry is retained independently of this condition status.",evidence.report.documentary_status));
                }
            }
            let receipt = verification::projection::project_digest_annotations(
                source,
                &geometry,
                &mut reference,
                &mut findings,
            )?;
            if let Some(capture) = capture {
                capture.record(
                    0,
                    receipt,
                    verification::ancestry::ProjectionKind::ReferenceCorrespondence,
                    vec![],
                    verification::ancestry::IntroductionKind::DeclaredLiteral,
                );
            }
            inherited.incomplete |= !inspection.scan_complete || !inspection.limits.is_empty();
            findings.push("Catalog cleavage uses actual tracked source strands and retains terminal groups, nicks and positional modifications. Complete cleavage is a declared reaction assumption, not measured enzyme completion or methylation permission; fragment choice is not proof of isolation.".into());
            (reference, fragment.duplex, findings)
        }
    };
    reference
        .metadata
        .insert(CONTEXT_KEY.into(), revision.into());
    reference.metadata.insert("source:material_derivation:scope".into(),"Target-free native processing recipe. The detached reference is not reusable material; retain exact source events and replay the native record. No independent target was assessed.".into());
    tracked.validate()?;
    reference.validate()?;
    let local_count = local.len();
    let total = inherited
        .findings_total
        .checked_add(local_count)
        .ok_or(Error::TooLarge)?;
    let mut findings = local;
    findings.extend(inherited.findings);
    let (findings, findings_total, findings_truncated) =
        limitations(findings, total, inherited.findings_truncated);
    Ok(FreshMaterial {
        reaction_conditions: inherited.reaction_conditions,
        reference,
        tracked,
        incomplete: inherited.incomplete,
        findings,
        findings_total,
        findings_truncated,
    })
}
pub fn prepare_install(w: &Workspace, request: &Request) -> Result<PreparedInstallation, Error> {
    let context_revision = revision(w, request)?;
    let fingerprint = run_lineage::workspace_fingerprint(w)?;
    let snapshots = workspace_snapshots(w, request)?;
    let fresh = fresh(&snapshots, request, &context_revision)?;
    let histories = DocumentHistory::restore_all(w)?;
    let parent_histories = request
        .sources
        .iter()
        .map(|s| {
            let history = histories
                .iter()
                .find(|h| h.id == s.event.document_id)
                .ok_or_else(|| bad("source history absent"))?;
            Ok((history, s.event.event_id.as_str()))
        })
        .collect::<Result<Vec<_>, Error>>()?;
    let prefix = match request.derivation {
        Derivation::ExactPcr { .. } => "pcr",
        Derivation::CatalogDigest { .. } => "material-digest",
        Derivation::MismatchPcr { .. } => "mismatch-pcr",
        Derivation::EndRepair { .. } => "material-repair",
        Derivation::OverlapExtensionPcr { .. } => "oe-pcr",
        Derivation::TerminalProcessing { .. } => "terminal-processing",
        Derivation::AnnealedOligos { .. } => "annealed-oligos",
    };
    let prefix = format!("{prefix}-{}", &fresh.reference.revision()[..24]);
    let mut id = prefix.clone();
    let mut suffix = 1usize;
    while w.documents.iter().any(|d| d.id == id) {
        id = format!("{prefix}-{suffix}");
        suffix = suffix.checked_add(1).ok_or(Error::TooLarge)?
    }
    let derived = DocumentHistory::derive_from_events(
        &parent_histories,
        id.clone(),
        fresh.reference.clone(),
        EVENT_OPERATION,
        AnalysisOrigin::context(
            AnalysisProducer::MaterialDerivation,
            context_revision.clone(),
        ),
    )?;
    let product = EventBinding {
        document_id: id.clone(),
        event_id: derived.head_event().into(),
        revision: derived.head_revision().into(),
    };
    let modifications = fresh.tracked.summary()?;
    let record = Record {
        schema_version: 5,
        engine_version: ENGINE_VERSION.into(),
        scientific_engine: verification::ENGINE_VERSION.into(),
        source_workspace_fingerprint: fingerprint,
        context_revision: context_revision.clone(),
        request: request.clone(),
        product: product.clone(),
        product_state_revision: modifications.state_revision.clone(),
        incomplete: fresh.incomplete,
    };
    let bytes = save_record(&record)?;
    let record_hash = blake3::hash(&bytes).to_hex().to_string();
    let mut next = w.clone();
    let primer_attachment = if request.save_primers {
        let table = declared_primer_table(&request.derivation)?
            .ok_or_else(|| bad("Only newly performed PCR has newly declared primers"))?;
        let table_bytes = encoded(&table, MAX_RECORD_BYTES)?;
        let hash = blake3::hash(&table_bytes).to_hex().to_string();
        if let Some(existing) = next.attachments.get(&hash) {
            if existing.kind != "oligo-table-v1"
                || existing.media_type != "application/json"
                || existing.bytes != table_bytes
            {
                return Err(bad("Primer table attachment collision."));
            }
        } else {
            if next.attachments.len() + 1 >= crate::documents::workspace::MAX_ATTACHMENTS {
                return Err(Error::TooLarge);
            }
            next.attachments.insert(
                hash.clone(),
                WorkspaceAttachment {
                    kind: "oligo-table-v1".into(),
                    name: format!("PCR primers: {}", request.name),
                    media_type: "application/json".into(),
                    bytes: table_bytes,
                },
            );
        }
        Some(hash)
    } else {
        None
    };
    derived.contribute(&mut next)?;
    next.active_document = id;
    if next
        .attachments
        .insert(
            record_hash.clone(),
            WorkspaceAttachment {
                kind: WORKSPACE_ATTACHMENT_KIND.into(),
                name: format!("Native material derivation: {}", request.name),
                media_type: MEDIA_TYPE.into(),
                bytes,
            },
        )
        .is_some()
    {
        return Err(bad("native derivation attachment collision"));
    }
    next.validate()?;
    let receipt = Receipt {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        record_hash,
        primer_attachment,
        sources: request.source_bindings(),
        product,
        incomplete: fresh.incomplete,
        findings: fresh.findings,
        findings_total: fresh.findings_total,
        findings_truncated: fresh.findings_truncated,
        duplex: fresh.tracked.duplex.summary()?,
        modifications,
        reaction_conditions: fresh.reaction_conditions,
    };
    receipt.validate_archive(&next)?;
    Ok(PreparedInstallation {
        workspace: next,
        receipt,
    })
}
pub fn install(w: &Workspace, request: &Request) -> Result<Installed, Error> {
    let prepared = prepare_install(w, request)?;
    Ok(Installed {
        archive: crate::documents::workspace::save_workspace(&prepared.workspace)?,
        receipt: prepared.receipt,
    })
}
fn record_for(w: &Workspace, product: &EventBinding) -> Result<(String, Record), Error> {
    run_material::event_shape(product)?;
    let mut found = None;
    let mut count = 0;
    for (hash, asset) in w
        .attachments
        .iter()
        .filter(|(_, a)| a.kind == WORKSPACE_ATTACHMENT_KIND)
    {
        count += 1;
        if count > MAX_RECORDS
            || asset.media_type != MEDIA_TYPE
            || blake3::hash(&asset.bytes).to_hex().as_str() != hash
        {
            return Err(bad("record count/media/hash differs"));
        }
        // Retained lookup must not restart ancestry validation from inside the
        // shared walker's attestation callback. The selected closed ancestry is
        // validated by that walker; public workspace operations validate w first.
        let record = load_record_shallow(&asset.bytes)?;
        if &record.product == product {
            if found.is_some() {
                return Err(bad("multiple records claim one material event"));
            }
            found = Some((hash.clone(), record));
        }
    }
    found.ok_or_else(|| bad("exact native derivation record unavailable"))
}
pub(crate) fn attest_record_shallow(
    w: &Workspace,
    record: &Record,
    source: Option<&[SourceSnapshot]>,
) -> Result<(), Error> {
    let (_, retained) = record_for(w, &record.product)?;
    if save_record_shallow(&retained)? != save_record_shallow(record)? {
        return Err(bad("closed record differs from exact retained bytes"));
    }
    if let Some(snapshots) = source {
        validate_source_snapshots(&record.request, snapshots)?;
    }
    for (i, bound) in record.request.sources.iter().enumerate() {
        let actual = run_material::resolve(w, &bound.event)?;
        if source.is_some_and(|snapshots| &snapshots[i].molecule != actual) {
            return Err(bad(
                "closed source snapshot differs from exact native event",
            ));
        }
        if matches!(record.request.derivation, Derivation::AnnealedOligos { .. }) {
            declared_oligo_source(actual)?;
        } else if matches!(bound.origin, SourceOrigin::OriginalReference) {
            original_reference_admission(actual)?;
        }
    }
    let product = run_material::resolve(w, &record.product)?;
    let event = w
        .history
        .iter()
        .find(|e| e.id == record.product.event_id)
        .ok_or_else(|| bad("product event absent"))?;
    let expected_origin = AnalysisOrigin::context(
        AnalysisProducer::MaterialDerivation,
        record.context_revision.clone(),
    );
    if event.transition
        != (Transition::Analysis {
            origin: expected_origin,
        })
        || event.operation != EVENT_OPERATION
        || event.parents != distinct_parent_events(&record.request)
        || product.metadata.get(CONTEXT_KEY) != Some(&record.context_revision)
        || !product.metadata.contains_key("duplex_model")
    {
        return Err(bad(
            "native material parent/operation/restrictive marker differs",
        ));
    }
    Ok(())
}
pub fn validate_workspace_records(w: &Workspace) -> Result<(), Error> {
    let mut products = BTreeSet::new();
    let mut count = 0;
    for (hash, asset) in w
        .attachments
        .iter()
        .filter(|(_, a)| a.kind == WORKSPACE_ATTACHMENT_KIND)
    {
        count += 1;
        if count > MAX_RECORDS
            || asset.media_type != MEDIA_TYPE
            || blake3::hash(&asset.bytes).to_hex().as_str() != hash
        {
            return Err(bad("record count/media/hash differs"));
        }
        let record = load_record(&asset.bytes)?;
        if !products.insert(record.product.event_id.clone()) {
            return Err(bad("duplicate derivation product owner"));
        }
        attest_record_shallow(w, &record, None)?;
        crate::material::capsule::validate_derived_request_origins(w, &record.request)?;
    }
    if w.history.iter().any(|event| {
        matches!(
            &event.transition,
            Transition::Analysis {
                origin: AnalysisOrigin {
                    producer: AnalysisProducer::MaterialDerivation,
                    record: crate::documents::history::transition::AnalysisRecord::Context { .. }
                }
            }
        ) && !products.contains(&event.id)
    }) {
        return Err(bad(
            "native material product event lacks its retained record",
        ));
    }
    Ok(())
}
impl Capsule {
    pub fn validate_documentary(&self) -> Result<(), Error> {
        crate::material::capsule::derived_capsule_preflight(self)
    }
}
/// Native record/source lookup after the caller's full Workspace validation;
/// no science, reserialization, source cloning or capsule copies.
pub(crate) fn documentary_capsule_parts(
    w: &Workspace,
    product: &EventBinding,
) -> Result<(String, Record), Error> {
    let (hash, record) = record_for(w, product)?;
    for source in &record.request.sources {
        collagene_model::bounds::model_preflight(&[run_material::resolve(w, &source.event)?])?;
    }
    run_material::resolve(w, product)?;
    Ok((hash, record))
}
pub fn capsule(w: &Workspace, product: &EventBinding) -> Result<Capsule, Error> {
    let mut roots = crate::material::capsule::documentary_origins(
        w,
        &[run_material::MaterialSelection::Derived(product.clone())],
        &[],
    )?
    .materials;
    let run_material::MaterialOrigin::Derived(result) = roots.remove(0) else {
        return Err(bad("material capsule kind differs"));
    };
    let _fresh = reconstruct(w, product)?;
    Ok(result)
}
pub fn validate_capsule_for_workspace(
    w: &Workspace,
    capsule: &Capsule,
    product: &EventBinding,
) -> Result<(), Error> {
    if &capsule.product != product {
        return Err(bad("capsule requires the exact saved material event"));
    }
    crate::material::capsule::validate_derived_capsule_origins(w, capsule)
}
pub fn replay(capsule: &Capsule, reference: &Molecule) -> Result<FreshMaterial, Error> {
    capsule.validate_documentary()?;
    let record = load_record_shallow(&capsule.record_bytes)?;
    replay_fields(&record, &capsule.source_snapshots, reference)
}
/// Capture actual current derivation projections and nested saved-source receipts.
pub(crate) fn replay_scoped(
    capsule: &Capsule,
    reference: &Molecule,
    replay: &mut verification::replay_budget::ReplayBudget<'_>,
    budget: Option<verification::ancestry::TraceBudget>,
) -> Result<
    (
        FreshMaterial,
        verification::ancestry::OriginBinding,
        Option<Result<verification::ancestry::ReconstructionTrace, Error>>,
    ),
    Error,
> {
    capsule.validate_documentary()?;
    let record = load_record_shallow(&capsule.record_bytes)?;
    collagene_model::bounds::model_preflight(&[reference])?;
    if record.scientific_engine != verification::ENGINE_VERSION {
        return Err(bad(
            "historical derivation requires current native reassessment",
        ));
    }
    let mut capture = budget.map(|budget| ancestry::Capture::new(budget, &record.context_revision));
    let fresh = fresh_observed(
        &capsule.source_snapshots,
        &record.request,
        &record.context_revision,
        &mut capture,
        replay,
    )?;
    validate_replayed(&record, reference, &fresh)?;
    let binding = verification::ancestry::OriginBinding::Material {
        product: record.product.clone(),
        record_hash: blake3::hash(&capsule.record_bytes).to_hex().to_string(),
        context_revision: record.context_revision.clone(),
        sources: record.request.source_bindings(),
    };
    let trace = capture.map(|capture| capture.finish(fresh.reference()));
    Ok((fresh, binding, trace))
}
fn validate_replayed(
    record: &Record,
    reference: &Molecule,
    fresh: &FreshMaterial,
) -> Result<(), Error> {
    if reference.revision() != record.product.revision
        || reference != fresh.reference()
        || fresh.incomplete != record.incomplete
        || fresh.tracked().summary()?.state_revision != record.product_state_revision
    {
        return Err(bad(
            "actual replay differs from saved material product/state/limitations",
        ));
    }
    Ok(())
}
fn replay_fields(
    record: &Record,
    source: &[SourceSnapshot],
    reference: &Molecule,
) -> Result<FreshMaterial, Error> {
    collagene_model::bounds::model_preflight(&[reference])?;
    if record.scientific_engine != verification::ENGINE_VERSION {
        return Err(bad(
            "historical record remains readable but requires current native reassessment before material reuse",
        ));
    }
    let fresh = fresh(source, &record.request, &record.context_revision)?;
    validate_replayed(record, reference, &fresh)?;
    Ok(fresh)
}
pub fn reconstruct(w: &Workspace, product: &EventBinding) -> Result<FreshMaterial, Error> {
    w.validate()?;
    let (_, record) = record_for(w, product)?;
    let source = workspace_snapshots(w, &record.request)?;
    let reference = run_material::resolve(w, product)?;
    replay_fields(&record, &source, reference)
}
fn summary_shape(summary: &ModificationSummary) -> Result<(), Error> {
    let total = summary
        .top_length
        .checked_add(summary.bottom_length)
        .ok_or(Error::TooLarge)?;
    let calls = summary
        .known_modified
        .checked_add(summary.known_unmodified)
        .and_then(|n| n.checked_add(summary.unknown_eligible))
        .ok_or(Error::TooLarge)?;
    if summary.engine_version != modified_duplex::ENGINE_VERSION
        || !hex(&summary.state_revision)
        || summary.top_length == 0
        || summary.bottom_length == 0
        || summary.top_length > collagene_model::MAX_SEQUENCE_LEN
        || summary.bottom_length > collagene_model::MAX_SEQUENCE_LEN
        || calls > total
        || summary.generation.len() > 1000
        || summary.generation.iter().any(|g| !hex(&g.source_revision))
    {
        return Err(bad("invalid bounded material summary"));
    }
    Ok(())
}

impl Receipt {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if self.engine_version != ENGINE_VERSION
            || !hex(expected)
            || self.context_revision != expected
            || !hex(&self.record_hash)
            || self
                .primer_attachment
                .as_ref()
                .is_some_and(|hash| !hex(hash))
            || self.findings.len() != self.findings_total.min(128)
            || self.findings.iter().any(|s| s.chars().count() > 2048)
            || self.findings_total > crate::material::lineage::MAX_RUN_BYTES
            || (self.findings_total > self.findings.len() && !self.findings_truncated)
            || self.duplex.top_length != self.modifications.top_length
            || self.duplex.bottom_length != self.modifications.bottom_length
        {
            return Err(bad("receipt context/bounds differ"));
        }
        reaction_conditions::validate(&self.reaction_conditions)?;
        validate_source_bindings(&self.sources)?;
        run_material::event_shape(&self.product)?;
        summary_shape(&self.modifications)?;
        encoded(self, MAX_RECORD_BYTES)?;
        Ok(())
    }
    pub fn validate_archive(&self, w: &Workspace) -> Result<(), Error> {
        self.validate_binding(&self.context_revision)?;
        w.validate()?;
        let (hash, record) = record_for(w, &self.product)?;
        if hash != self.record_hash
            || record.context_revision != self.context_revision
            || record.request.source_bindings() != self.sources
            || record.incomplete != self.incomplete
            || record.product_state_revision != self.modifications.state_revision
        {
            return Err(bad("receipt differs from exact installed record"));
        }
        // Documentary association only: current scientific claims use fresh
        // prepare/reconstruct, never this imported receipt or its retained hashes.
        let source_bindings = record.source_bindings()?;
        if self.reaction_conditions.iter().any(|evidence| {
            !source_bindings
                .iter()
                .any(|binding| binding.revision == evidence.source_revision)
        }) || (!self.incomplete
            && self
                .reaction_conditions
                .iter()
                .any(|evidence| !evidence.satisfies_request()))
        {
            return Err(bad(
                "receipt documentary evidence is not associated with retained source ancestry",
            ));
        }
        if let Derivation::CatalogDigest {
            selection,
            conditions: Some(options),
            ..
        } = &record.request.derivation
        {
            let source = &record.request.template()?.event.revision;
            let evidence = self
                .reaction_conditions
                .iter()
                .find(|evidence| {
                    &evidence.source_revision == source
                        && evidence.report.selection == *selection
                        && evidence.report.options == *options
                })
                .ok_or_else(|| bad("receipt drops the declared documentary reaction setup"))?;
            evidence.report.validate_for(
                &catalog_digest::builtin_catalog()?,
                selection,
                options,
            )?;
        }
        if self.primer_attachment.is_some() != record.request.save_primers {
            return Err(bad(
                "Installation primer-table choice differs from its bound request.",
            ));
        }
        if let Some(hash) = &self.primer_attachment {
            let table = declared_primer_table(&record.request.derivation)?
                .ok_or_else(|| bad("A non-PCR receipt cannot attach new primers"))?;
            let expected = encoded(&table, MAX_RECORD_BYTES)?;
            let asset = w
                .attachments
                .get(hash)
                .ok_or_else(|| bad("Installation primer table unavailable."))?;
            if asset.kind != "oligo-table-v1"
                || asset.media_type != "application/json"
                || asset.bytes != expected
                || blake3::hash(&asset.bytes).to_hex().as_str() != hash
            {
                return Err(bad(
                    "Installation primer table differs from the full declared new oligos.",
                ));
            }
        }
        Ok(())
    }
}
impl Installed {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.receipt.validate_binding(expected)?;
        let w = crate::documents::workspace::load_workspace(&self.archive)?;
        self.receipt.validate_archive(&w)
    }
}
/// Literal primers remain recipe data even when a user chooses no order export.
pub fn primer_table(record: &Record) -> Result<Option<crate::formats::oligos::OligoTable>, Error> {
    record.validate_documentary()?;
    declared_primer_table(&record.request.derivation)
}
fn declared_primer_table(
    derivation: &Derivation,
) -> Result<Option<crate::formats::oligos::OligoTable>, Error> {
    match derivation {
        Derivation::ExactPcr {
            forward, reverse, ..
        }
        | Derivation::MismatchPcr {
            forward, reverse, ..
        } => crate::cloning::pcr::outputs::primer_table(forward, reverse).map(Some),
        Derivation::OverlapExtensionPcr { preparation, .. } => {
            preparation.validate_shape()?;
            let primers = [&preparation.left, &preparation.right, &preparation.outer]
                .into_iter()
                .flat_map(|pair| [pair.forward.clone(), pair.reverse.clone()])
                .collect::<Vec<_>>();
            crate::planning::orders::from_primers(&primers).map(Some)
        }
        _ => Ok(None),
    }
}
fn source_origin_preflight(source: &Molecule, origin: &SourceOrigin) -> Result<(), Error> {
    collagene_model::bounds::model_preflight(&[source])?;
    source.validate()?;
    match origin {
        SourceOrigin::OriginalReference => {
            original_source(source)?;
        }
        SourceOrigin::Saved(origin) => {
            origin.validate_documentary()?;
            if origin.event().revision != source.revision() {
                return Err(bad("preview source differs from exact saved material"));
            }
        }
    }
    Ok(())
}
fn preview_receipt_shape(
    engine: &str,
    context: &str,
    expected: &str,
    source: &str,
    findings: &[String],
    total: usize,
    truncated: bool,
) -> Result<(), Error> {
    if engine != ENGINE_VERSION
        || context != expected
        || !hex(expected)
        || !hex(source)
        || findings.len() != total.min(128)
        || total > crate::material::lineage::MAX_RUN_BYTES
        || findings.iter().any(|s| s.chars().count() > 2048)
        || (total > findings.len() && !truncated)
    {
        return Err(bad("processing preview binding/excerpt differs"));
    }
    Ok(())
}
fn preview_source_shape(
    duplex: &chemistry::DuplexSummary,
    modifications: &ModificationSummary,
) -> Result<(), Error> {
    summary_shape(modifications)?;
    if duplex.top_length != modifications.top_length
        || duplex.bottom_length != modifications.bottom_length
        || duplex.model.len() > 256
        || duplex.nicks.len() > 10000
        || duplex.bottom_start.unsigned_abs() > collagene_model::MAX_SEQUENCE_LEN as u64
        || duplex.ends.as_ref().is_some_and(|(a, b)| {
            a.overhang.len() > collagene_model::MAX_SEQUENCE_LEN
                || b.overhang.len() > collagene_model::MAX_SEQUENCE_LEN
        })
    {
        return Err(bad(
            "processing preview source geometry/modifications differ",
        ));
    }
    Ok(())
}
fn workspace_source<'a>(
    w: &'a Workspace,
    source: &EventBinding,
    origin: &SourceOrigin,
) -> Result<&'a Molecule, Error> {
    origin.validate_shallow_for(source)?;
    crate::material::capsule::material_origins_preflight(match origin {
        SourceOrigin::OriginalReference => None,
        SourceOrigin::Saved(origin) => Some(origin),
    })?;
    w.validate()?;
    if let SourceOrigin::Saved(origin) = origin {
        origin.validate_for_workspace(w, source)?;
    }
    let molecule = run_material::resolve(w, source)?;
    source_origin_preflight(molecule, origin)?;
    Ok(molecule)
}

fn validate_source_bindings(sources: &[SourceBinding]) -> Result<(), Error> {
    let expected: &[SourceRole] = match sources.len() {
        1 => &[SourceRole::Template],
        2 => &[SourceRole::TopOligo, SourceRole::BottomOligo],
        _ => return Err(bad("unsupported source count")),
    };
    for (source, role) in sources.iter().zip(expected) {
        if source.role != *role {
            return Err(bad("source roles changed canonical order"));
        }
        run_material::event_shape(&source.event)?;
    }
    Ok(())
}
fn distinct_parent_events(request: &Request) -> Vec<String> {
    let mut parents = Vec::new();
    for source in &request.sources {
        if !parents.contains(&source.event.event_id) {
            parents.push(source.event.event_id.clone())
        }
    }
    parents
}
fn declared_oligo_source(molecule: &Molecule) -> Result<(), Error> {
    collagene_model::bounds::model_preflight(&[molecule])?;
    verification::terminal::declared_oligo(molecule)?;
    if molecule.sequence.is_empty() || molecule.sequence.len() > chemistry::MAX_ANNEAL_OLIGO_BASES {
        return Err(Error::TooLarge);
    }
    if !molecule
        .sequence
        .as_str()
        .bytes()
        .all(|b| b"ACGT".contains(&b))
    {
        return Err(bad("standard annealing oligos require concrete DNA"));
    }
    Ok(())
}
pub(crate) fn validate_source_snapshots(
    request: &Request,
    snapshots: &[SourceSnapshot],
) -> Result<(), Error> {
    request.validate_shallow()?;
    if snapshots.len() != request.sources.len() {
        return Err(bad("source snapshot count differs from retained roles"));
    }
    for (snapshot, source) in snapshots.iter().zip(&request.sources) {
        if snapshot.role != source.role || snapshot.molecule.revision() != source.event.revision {
            return Err(bad(
                "source snapshot role/order or revision differs from exact source event",
            ));
        }
        collagene_model::bounds::model_preflight(&[&snapshot.molecule])?;
        snapshot.molecule.validate()?;
        if matches!(request.derivation, Derivation::AnnealedOligos { .. }) {
            declared_oligo_source(&snapshot.molecule)?;
        } else if matches!(source.origin, SourceOrigin::OriginalReference) {
            original_reference_admission(&snapshot.molecule)?;
        }
    }
    Ok(())
}
fn workspace_snapshots(w: &Workspace, request: &Request) -> Result<Vec<SourceSnapshot>, Error> {
    request
        .sources
        .iter()
        .map(|source| {
            Ok(SourceSnapshot {
                role: source.role,
                molecule: run_material::resolve(w, &source.event)?.clone(),
            })
        })
        .collect()
}
fn fresh(
    snapshots: &[SourceSnapshot],
    request: &Request,
    revision: &str,
) -> Result<FreshMaterial, Error> {
    let mut cleavage = modified_duplex::CleavageBudget::new(100_000_000)?;
    let mut replay = verification::replay_budget::ReplayBudget::plan_local(&mut cleavage);
    fresh_observed(snapshots, request, revision, &mut None, &mut replay)
}
fn fresh_observed(
    snapshots: &[SourceSnapshot],
    request: &Request,
    revision: &str,
    capture: &mut Option<ancestry::Capture>,
    replay: &mut verification::replay_budget::ReplayBudget<'_>,
) -> Result<FreshMaterial, Error> {
    validate_source_snapshots(request, snapshots)?;
    if let Derivation::AnnealedOligos { options, selection } = &request.derivation {
        let mut findings = Vec::new();
        let preview = verification::terminal::simulate_oligo_annealing(
            &snapshots[0].molecule,
            &snapshots[1].molecule,
            options,
        )?;
        if annealing::anneal_context(&request.source_bindings(), options)?
            != selection.expected_preview_revision
        {
            return Err(bad("stale exact-source annealing preview selection"));
        }
        let index = preview
            .annealing
            .candidates
            .iter()
            .position(|c| c.id == selection.candidate_id)
            .ok_or_else(|| bad("selected native annealing candidate is unavailable"))?;
        let projected = verification::terminal::annealed_document_inner(
            &snapshots[0].molecule,
            &snapshots[1].molecule,
            &preview.annealing,
            index,
            &mut findings,
        )?;
        let mut reference = projected.molecule;
        if let Some(capture) = capture {
            for snapshot in snapshots {
                capture.root(snapshot.role, &snapshot.molecule, None)?;
            }
            for (input, receipt) in projected.projections.into_iter().enumerate() {
                capture.record(
                    input,
                    receipt,
                    verification::ancestry::ProjectionKind::ReferenceCorrespondence,
                    vec![],
                    verification::ancestry::IntroductionKind::DeclaredLiteral,
                );
            }
        }
        annealing::protect_anneal_preview(&mut reference)?;
        if reference.revision() != selection.expected_product_revision {
            return Err(bad("selected native annealed product changed"));
        }
        let tracked =
            TrackedDuplex::annealed_unmodified(&preview.annealing.candidates[index].duplex)?;
        let incomplete = !preview.annealing.complete || !preview.annealing.unique;
        reference.name = request.name.clone();
        reference
            .metadata
            .insert(CONTEXT_KEY.into(), revision.into());
        reference.metadata.insert("source:material_derivation:scope".into(),"Target-free native standard-oligo annealing. Retain both exact role-bound sources and recipe; no intended target or experimental isolation was assessed.".into());
        findings.extend(preview.annealing.warnings);
        findings.push("Both inputs are explicitly declared independently supplied, standard unmodified ssDNA oligos. End groups are the recorded chemistry hypotheses. No saved duplex was denatured, and no primer order, experimental isolation or intended target was inferred.".into());
        if incomplete {
            findings.push("Selected annealed material remains conditional: alternative geometries or incomplete offset enumeration remain unresolved by candidate selection.".into());
        }
        let total = findings.len();
        let (findings, findings_total, findings_truncated) = limitations(findings, total, false);
        reference.validate()?;
        tracked.validate()?;
        Ok(FreshMaterial {
            reaction_conditions: vec![],
            reference,
            tracked,
            incomplete,
            findings,
            findings_total,
            findings_truncated,
        })
    } else {
        fresh_template(&snapshots[0].molecule, request, revision, capture, replay)
    }
}
