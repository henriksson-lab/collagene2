//! Target-free native PCR, reciprocal extension, catalog-digest and end-repair derivations. Opaque, bounded recipe
//! packages break typed recursion; only fresh native execution yields material.
use crate::{
    Alphabet, Error, Molecule,
    analysis::{self, Primer},
    catalog_digest::{self, CatalogSelection},
    chemistry,
    history::DocumentHistory,
    modified_duplex::{self, ModificationSummary, TrackedDuplex},
    restriction_inspection::{self, InspectionOptions, RestrictionInspection},
    run_lineage::{self, EventBinding},
    run_material::{self, MaterialOrigin},
    verification,
    workspace::{Workspace, WorkspaceAttachment},
};
use base64::{Engine as _, engine::general_purpose::STANDARD};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};
pub const ENGINE_VERSION: &str = "material-derivation-2";
pub const WORKSPACE_ATTACHMENT_KIND: &str = "material-derivation-v2";
pub const MEDIA_TYPE: &str = "application/vnd.collagene.material-derivation";
pub const EVENT_OPERATION: &str = "Target-free native material derivation v2";
pub const MAX_RECORD_BYTES: usize = 8 * 1024 * 1024;
pub const MAX_CAPSULE_BYTES: usize = 32 * 1024 * 1024;
pub const MAX_RECORDS: usize = 128;
pub const MAX_SOURCES: usize = 2;
const MAGIC: &[u8; 8] = b"MDRV\0\0\0\x02";
const HEADER: usize = 44;
const CONTEXT_KEY: &str = "source:material_derivation:context";
macro_rules! dto{($($t:tt)*)=>{#[derive(Debug,Clone,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {#[derive(Copy,PartialEq,Eq,PartialOrd,Ord)]#[serde(rename_all="snake_case")]pub enum SourceRole{Template,TopOligo,BottomOligo}}
dto! {pub struct BoundSource{pub role:SourceRole,pub event:EventBinding,pub origin:SourceOrigin}}
dto! {#[derive(PartialEq,Eq)]pub struct SourceBinding{pub role:SourceRole,pub event:EventBinding}}
dto! {pub struct SourceSnapshot{pub role:SourceRole,pub molecule:Molecule}}
dto! {pub enum SourceOrigin{OriginalReference,Saved(MaterialOrigin)}}
dto! {pub struct PcrSelection{pub selected_index:usize,pub expected_preview_revision:String,#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]pub geometry:analysis::PcrProduct}}
dto! {pub struct MismatchPcrSelection{pub selected_index:usize,pub expected_preview_revision:String,#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]pub geometry:crate::mutagenesis::MutantPcrProduct}}
dto! {pub struct OverlapExtensionSelection{pub expected_preview_revision:String,pub expected_product_revision:String}}
dto! {#[derive(Copy,PartialEq,Eq)]pub enum TerminalAction{TreatEnds{treatment:chemistry::EndTreatment},Tail3Prime{base:crate::tailing::TailingBase}}}
dto! {pub struct AnnealSelection{pub expected_preview_revision:String,pub candidate_id:String,pub expected_product_revision:String}}
dto! {pub enum Derivation{ExactPcr{forward:Primer,reverse:Primer,template_strand:chemistry::Strand,selection:PcrSelection},CatalogDigest{selection:CatalogSelection,options:InspectionOptions,fragment_index:usize},MismatchPcr{forward:Primer,reverse:Primer,options:crate::mutagenesis::PcrOptions,selection:MismatchPcrSelection},EndRepair{method:chemistry::EndRepairMethod},OverlapExtensionPcr{preparation:crate::overlap_extension_pcr::Request,selection:OverlapExtensionSelection},TerminalProcessing{action:TerminalAction},AnnealedOligos{options:chemistry::AnnealOptions,selection:AnnealSelection}}}
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
dto! {pub struct Receipt{pub engine_version:String,pub context_revision:String,pub record_hash:String,pub primer_attachment:Option<String>,pub sources:Vec<SourceBinding>,pub product:EventBinding,pub incomplete:bool,pub findings:Vec<String>,pub findings_total:usize,pub findings_truncated:bool,pub duplex:chemistry::DuplexSummary,pub modifications:ModificationSummary}}
pub struct PreparedInstallation {
    pub workspace: Workspace,
    pub receipt: Receipt,
}
dto! {pub struct Installed{pub archive:Vec<u8>,pub receipt:Receipt}}
dto! {pub struct DigestPreview{pub engine_version:String,pub context_revision:String,pub selection:CatalogSelection,pub inspection:RestrictionInspection,pub source_duplex:chemistry::DuplexSummary,pub source_modifications:ModificationSummary,pub fragment_modifications:Vec<ModificationSummary>,pub incomplete:bool,pub findings:Vec<String>,pub findings_total:usize,pub findings_truncated:bool}}
/// Not serializable, not publicly constructible, and never reconstructed from summaries.
#[derive(Debug)]
pub struct FreshMaterial {
    reference: Molecule,
    tracked: TrackedDuplex,
    incomplete: bool,
    findings: Vec<String>,
    findings_total: usize,
    findings_truncated: bool,
}
impl FreshMaterial {
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
                anneal_options_shape(options)?;
                if !hex(&selection.expected_preview_revision)
                    || !hex(&selection.expected_product_revision)
                    || !candidate_id_shape(&selection.candidate_id)
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
                    crate::pcr_outputs::primer(p)?;
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
                    || selection.geometry.sequence.len() > crate::model::MAX_SEQUENCE_LEN
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
                crate::assembly_recipe::mismatch_pcr_shape(forward, reverse, options)?;
                if selection.selected_index >= 10000
                    || !hex(&selection.expected_preview_revision)
                    || selection.geometry.sequence.is_empty()
                    || selection.geometry.sequence.len() > crate::model::MAX_SEQUENCE_LEN
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
                        != overlap_extension_context(
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
                fragment_index,
            } => {
                if self.save_primers {
                    return Err(bad(
                        "Catalog digest cannot create a new PCR primer order; existing upstream primers remain source history.",
                    ));
                }
                selection.validate()?;
                options.validate()?;
                if *fragment_index >= options.max_fragments {
                    return Err(bad("fragment index exceeds declared bounds"));
                }
            }
        }
        Ok(())
    }
}
fn context(fingerprint: &str, request: &Request, science: &str) -> Result<String, Error> {
    crate::assembly_recipe::hash(&(ENGINE_VERSION, science, fingerprint, request))
}
impl Record {
    pub fn validate_documentary(&self) -> Result<(), Error> {
        crate::run_material_capsule::derived_record_preflight(self)
    }
    pub(crate) fn validate_shallow(&self) -> Result<(), Error> {
        if self.schema_version != 2
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
pub(crate) fn load_record_shallow(bytes: &[u8]) -> Result<Record, Error> {
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
    crate::model::bounds::model_preflight(&[source])?;
    source.validate()?;
    crate::format_fidelity::reference_export_guard(source)?;
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
    crate::model::bounds::model_preflight(&[source])?;
    source.validate()?;
    let (tracked,incomplete,findings,total,truncated)=match origin {
        SourceOrigin::OriginalReference=>(original_source(source)?,false,vec!["Original source is an explicit intact DNA reference hypothesis with unknown source modifications and terminal groups; no measured stock chemistry or amplification success is inferred.".into()],1,false),
        SourceOrigin::Saved(origin)=>{let (_,tracked,incomplete,findings,total,truncated)=origin.replay(source)?.into_parts();(tracked,incomplete,findings,total,truncated)},
    };
    let (findings, findings_total, findings_truncated) = limitations(findings, total, truncated);
    Ok(FreshMaterial {
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
    crate::model::bounds::model_preflight(&[source])?;
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
        if crate::assembly_recipe::mismatch_pcr_work(
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
        crate::overlap_extension_execution::preflight(source, preparation)?;
    }
    if let SourceOrigin::OriginalReference = &template.origin {
        original_source(source)?;
    }
    Ok(())
}
pub fn preflight(w: &Workspace, request: &Request) -> Result<(), Error> {
    request.validate_shallow()?;
    w.validate()?;
    crate::run_material_capsule::validate_derived_request_origins(w, request)?;
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
    if w.attachments.len() >= crate::workspace::MAX_ATTACHMENTS
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
) -> Result<FreshMaterial, Error> {
    let template = request.template()?;
    request_source(source, request)?;
    let mut inherited = fresh_source(source, &template.origin)?;
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
            let previews = match &template.origin {
                SourceOrigin::OriginalReference => verification::simulate_pcr_documents(
                    source,
                    forward,
                    reverse,
                    *template_strand,
                )?,
                SourceOrigin::Saved(origin) => crate::material_pcr::simulate_from_fresh(
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
                )?,
            };
            if previews.len() != 1 || selection.selected_index != 0 {
                return Err(bad(
                    "initial reusable PCR profile requires exactly one complete geometry; multiple products remain reference previews, never inferred isolation",
                ));
            }
            let preview = &previews[0];
            if preview.molecule.revision() != selection.expected_preview_revision
                || crate::assembly_recipe::hash(&preview.geometry)?
                    != crate::assembly_recipe::hash(&selection.geometry)?
            {
                return Err(bad("selected PCR preview/geometry changed"));
            }
            let frame =
                crate::pcr_source::Template::fresh(source, &inherited.tracked, *template_strand)?;
            let geometries = frame.geometries(forward, reverse)?;
            if geometries.len() != 1 {
                return Err(bad("actual exact PCR geometry is not unique"));
            }
            let geometry = &geometries[0];
            let tracked = TrackedDuplex::exact_pcr(&geometry.product.sequence)?;
            let mut reference = tracked.duplex.molecule(&request.name)?;
            let mut findings = Vec::new();
            verification::project_pcr_frame_annotations(
                frame,
                geometry,
                &mut reference,
                &mut findings,
            )?;
            findings.push("Unique exact PCR produces newly synthesized unmodified blunt duplex strands with 5′OH/3′OH in the existing ideal-PCR model. Primer chemistry, yields, cleanup, amplification efficiency and experimental isolation are not measured.".into());
            (reference, tracked, findings)
        }
        Derivation::MismatchPcr {
            forward,
            reverse,
            options,
            selection,
        } => {
            let preview = mismatch_from_fresh(
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
                || crate::assembly_recipe::hash(geometry)?
                    != crate::assembly_recipe::hash(&selection.geometry)?
            {
                return Err(bad("selected mismatch-PCR preview/geometry changed"));
            }
            inherited.incomplete |= !preview.preview.report.enumeration_complete;
            let tracked = TrackedDuplex::exact_pcr(&geometry.sequence)?;
            let mut reference = tracked.duplex.molecule(&request.name)?;
            let mut findings = preview.preview.report.warnings.clone();
            verification::project_mutant_pcr(source, geometry, &mut reference, &mut findings)?;
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
            let executed = crate::overlap_extension_execution::execute(
                source,
                Some(&inherited.tracked),
                preparation,
            )?;
            let mut product=executed.product.ok_or_else(||bad("Reusable OE material requires one actual geometry in each source PCR and outer PCR; no product is inferred from a desired sequence"))?;
            if product.revision() != selection.expected_product_revision {
                return Err(bad(
                    "Selected outer-PCR product changed; inspect the source-backed OE preparation again",
                ));
            }
            product.name = request.name.clone();
            inherited.incomplete |= executed.incomplete;
            let tracked = TrackedDuplex::exact_pcr(&product.sequence)?;
            let mut findings = executed.findings;
            findings.push("Three literal primer pairs produce the actual outer-PCR product. Newly synthesized unmodified blunt strands carry 5′OH/3′OH under the ideal-PCR model. Source/terminal uncertainty remains inherited; amplification efficiency and experimental isolation are not established, and no intended target was checked.".into());
            (product, tracked, findings)
        }
        Derivation::TerminalProcessing { action } => {
            terminal_from_fresh(source, &inherited, action, &request.name)?
        }
        Derivation::EndRepair { method } => {
            let repair = chemistry::repair_ends(&inherited.tracked.duplex, *method)?;
            let tracked = inherited.tracked.repair(*method)?;
            if tracked.duplex != repair.duplex {
                return Err(bad("native repair geometry and tracked mappings differ"));
            }
            let mut findings = repair.warnings.clone();
            let reference = verification::repaired_document(
                source,
                &repair,
                request.name.clone(),
                &mut findings,
            )?;
            findings.push("Direct end repair uses actual source strands and template-supported fill/removal only. Surviving positional modifications, original 5′ groups and internal nicks remain; no phosphorylation, nick sealing, supplier-vector chemistry or target success is inferred.".into());
            (reference, tracked, findings)
        }
        Derivation::CatalogDigest {
            selection,
            options,
            fragment_index,
        } => {
            let catalog = catalog_digest::builtin_catalog()?;
            selection.validate_for(&catalog)?;
            let enzymes = selection.resolved_enzymes(&catalog)?;
            let role = crate::restriction_history::source_role(&template.event)?;
            let mut selected = None;
            let inspection = restriction_inspection::inspect_with_digest_for_source(
                &role,
                source,
                &enzymes,
                options,
                &inherited.tracked.duplex.summary()?,
                || {
                    let mut result = Vec::new();
                    for (index, fragment) in modified_duplex::digest(&inherited.tracked, &enzymes)?
                        .into_iter()
                        .enumerate()
                    {
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
            verification::project_digest_annotations(
                source,
                &geometry,
                &mut reference,
                &mut findings,
            )?;
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
    )?;
    let product = EventBinding {
        document_id: id.clone(),
        event_id: derived.head_event().into(),
        revision: derived.head_revision().into(),
    };
    let modifications = modified_duplex::summary(&fresh.tracked)?;
    let record = Record {
        schema_version: 2,
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
            if next.attachments.len() + 1 >= crate::workspace::MAX_ATTACHMENTS {
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
        archive: crate::workspace::save_workspace(&prepared.workspace)?,
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
    if event.operation != EVENT_OPERATION
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
        crate::run_material_capsule::validate_derived_request_origins(w, &record.request)?;
    }
    if w.history
        .iter()
        .any(|e| e.operation == EVENT_OPERATION && !products.contains(&e.id))
    {
        return Err(bad(
            "native material product event lacks its retained record",
        ));
    }
    Ok(())
}
impl Capsule {
    pub fn validate_documentary(&self) -> Result<(), Error> {
        crate::run_material_capsule::derived_capsule_preflight(self)
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
        crate::model::bounds::model_preflight(&[run_material::resolve(w, &source.event)?])?;
    }
    run_material::resolve(w, product)?;
    Ok((hash, record))
}
pub fn capsule(w: &Workspace, product: &EventBinding) -> Result<Capsule, Error> {
    let mut roots = run_material::documentary_origins(
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
    crate::run_material_capsule::validate_derived_capsule_origins(w, capsule)
}
pub fn replay(capsule: &Capsule, reference: &Molecule) -> Result<FreshMaterial, Error> {
    capsule.validate_documentary()?;
    let record = load_record_shallow(&capsule.record_bytes)?;
    replay_fields(&record, &capsule.source_snapshots, reference)
}
fn replay_fields(
    record: &Record,
    source: &[SourceSnapshot],
    reference: &Molecule,
) -> Result<FreshMaterial, Error> {
    crate::model::bounds::model_preflight(&[reference])?;
    if record.scientific_engine != verification::ENGINE_VERSION {
        return Err(bad(
            "historical record remains readable but requires current native reassessment before material reuse",
        ));
    }
    let fresh = fresh(source, &record.request, &record.context_revision)?;
    if reference.revision() != record.product.revision
        || reference != fresh.reference()
        || fresh.incomplete != record.incomplete
        || modified_duplex::summary(fresh.tracked())?.state_revision
            != record.product_state_revision
    {
        return Err(bad(
            "actual replay differs from saved material product/state/limitations",
        ));
    }
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
        || summary.top_length > crate::model::MAX_SEQUENCE_LEN
        || summary.bottom_length > crate::model::MAX_SEQUENCE_LEN
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
            || self.findings_total > crate::run_lineage::MAX_RUN_BYTES
            || (self.findings_total > self.findings.len() && !self.findings_truncated)
            || self.duplex.top_length != self.modifications.top_length
            || self.duplex.bottom_length != self.modifications.bottom_length
        {
            return Err(bad("receipt context/bounds differ"));
        }
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
        let w = crate::workspace::load_workspace(&self.archive)?;
        self.receipt.validate_archive(&w)
    }
}
/// Literal primers remain recipe data even when a user chooses no order export.
pub fn primer_table(record: &Record) -> Result<Option<crate::oligo_formats::OligoTable>, Error> {
    record.validate_documentary()?;
    declared_primer_table(&record.request.derivation)
}
fn declared_primer_table(
    derivation: &Derivation,
) -> Result<Option<crate::oligo_formats::OligoTable>, Error> {
    match derivation {
        Derivation::ExactPcr {
            forward, reverse, ..
        }
        | Derivation::MismatchPcr {
            forward, reverse, ..
        } => crate::pcr_outputs::primer_table(forward, reverse).map(Some),
        Derivation::OverlapExtensionPcr { preparation, .. } => {
            preparation.validate_shape()?;
            let primers = [&preparation.left, &preparation.right, &preparation.outer]
                .into_iter()
                .flat_map(|pair| [pair.forward.clone(), pair.reverse.clone()])
                .collect::<Vec<_>>();
            crate::oligo_orders::from_primers(&primers).map(Some)
        }
        _ => Ok(None),
    }
}

pub fn source_bindings(record: &Record) -> Result<Vec<EventBinding>, Error> {
    record.validate_documentary()?;
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
    for source in &record.request.sources {
        add(&source.event);
    }
    // Collect callback results separately so only one closure mutably borrows the map.
    let mut runs = Vec::new();
    let mut derived = Vec::new();
    let mut graphs = Vec::new();
    crate::run_material_capsule::derived_ancestry_walk(
        record,
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
            if r.product != record.product {
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
pub fn required_run_ids(record: &Record) -> Result<Vec<String>, Error> {
    let mut ids = BTreeSet::new();
    crate::run_material_capsule::derived_ancestry_walk(
        record,
        |origin, _| {
            ids.insert(origin.selection.run_id.clone());
            Ok(())
        },
        |_, _| Ok(()),
        |_, _| Ok(()),
    )?;
    Ok(ids.into_iter().collect())
}
pub fn required_derivation_products(record: &Record) -> Result<Vec<EventBinding>, Error> {
    let mut ids = BTreeMap::new();
    crate::run_material_capsule::derived_ancestry_walk(
        record,
        |_, _| Ok(()),
        |r, _| {
            if r.product != record.product {
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
    record: &Record,
) -> Result<Vec<crate::topo_graph_history::Selection>, Error> {
    let mut records = BTreeMap::new();
    crate::run_material_capsule::derived_ancestry_walk(
        record,
        |_, _| Ok(()),
        |_, _| Ok(()),
        |graph, _| {
            let bytes = crate::topo_graph_history::save_record(graph)?;
            let hash = blake3::hash(&bytes).to_hex().to_string();
            records.insert(
                hash.clone(),
                crate::topo_graph_history::Selection {
                    record_hash: hash,
                    state_revision: graph.state_revision.clone(),
                },
            );
            Ok(())
        },
    )?;
    Ok(records.into_values().collect())
}
/// Source-bound reference previews only; this helper creates no native product.
pub fn pcr_preview_revision(
    source: &Molecule,
    origin: &SourceOrigin,
    forward: &Primer,
    reverse: &Primer,
    template_strand: chemistry::Strand,
) -> Result<String, Error> {
    crate::model::bounds::model_preflight(&[source])?;
    source.validate()?;
    for primer in [forward, reverse] {
        if primer.annealing.len() > 1000 || primer.tail.len() > 1000 {
            return Err(Error::TooLarge);
        }
        crate::pcr_outputs::primer(primer)?;
    }
    match origin {
        SourceOrigin::OriginalReference => {
            if template_strand != chemistry::Strand::Top {
                return Err(bad(
                    "Original reference PCR uses explicit Top; Bottom requires native source replay",
                ));
            }
            original_source(source)?;
        }
        SourceOrigin::Saved(origin) => {
            origin.validate_documentary()?;
            if origin.event().revision != source.revision() {
                return Err(bad("PCR preview source differs from saved origin."));
            }
        }
    }
    crate::assembly_recipe::hash(&(
        ENGINE_VERSION,
        "pcr-preview",
        verification::ENGINE_VERSION,
        source,
        origin,
        forward,
        reverse,
        template_strand,
    ))
}
pub fn preview_pcr(
    source: &Molecule,
    origin: &SourceOrigin,
    forward: &Primer,
    reverse: &Primer,
    template_strand: chemistry::Strand,
) -> Result<Vec<verification::PcrDocument>, Error> {
    pcr_preview_revision(source, origin, forward, reverse, template_strand)?;
    match origin {
        SourceOrigin::OriginalReference => {
            verification::simulate_pcr_documents(source, forward, reverse, template_strand)
        }
        SourceOrigin::Saved(origin) => {
            crate::material_pcr::simulate(source, origin, forward, reverse, template_strand)
        }
    }
}
pub fn digest_inspection_revision(
    role: &str,
    source: &Molecule,
    origin: &SourceOrigin,
    selection: &CatalogSelection,
    options: &InspectionOptions,
) -> Result<String, Error> {
    crate::model::bounds::model_preflight(&[source])?;
    source.validate()?;
    options.validate()?;
    let catalog = catalog_digest::builtin_catalog()?;
    catalog_digest::source_preflight(&catalog, selection, role, source, options)?;
    match origin {
        SourceOrigin::OriginalReference => {
            original_source(source)?;
        }
        SourceOrigin::Saved(origin) => {
            origin.validate_documentary()?;
            if origin.event().revision != source.revision() {
                return Err(bad("inspection origin/source revision differs"));
            }
        }
    }
    crate::assembly_recipe::hash(&(
        ENGINE_VERSION,
        "digest-inspection",
        verification::ENGINE_VERSION,
        role,
        source,
        origin,
        selection,
        options,
    ))
}
pub fn inspect_digest(
    role: &str,
    source: &Molecule,
    origin: &SourceOrigin,
    selection: &CatalogSelection,
    options: &InspectionOptions,
) -> Result<DigestPreview, Error> {
    let context_revision = digest_inspection_revision(role, source, origin, selection, options)?;
    let fresh = fresh_source(source, origin)?;
    let catalog = catalog_digest::builtin_catalog()?;
    let enzymes = selection.resolved_enzymes(&catalog)?;
    let mut fragment_modifications = Vec::new();
    let inspection = restriction_inspection::inspect_with_digest_for_source(
        role,
        source,
        &enzymes,
        options,
        &fresh.tracked.duplex.summary()?,
        || {
            let mut result = Vec::new();
            for fragment in modified_duplex::digest(&fresh.tracked, &enzymes)? {
                fragment_modifications.push(modified_duplex::summary(&fragment.duplex)?);
                result.push(chemistry::DigestFragment {
                    duplex: fragment.duplex.duplex,
                    top_start: fragment.top_start,
                    bottom_start: fragment.bottom_start,
                })
            }
            Ok(result)
        },
    )?;
    if !inspection.products_determined {
        fragment_modifications.clear()
    }
    let report = DigestPreview {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        selection: selection.clone(),
        source_duplex: fresh.tracked.duplex.summary()?,
        source_modifications: modified_duplex::summary(&fresh.tracked)?,
        fragment_modifications,
        incomplete: fresh.incomplete || !inspection.scan_complete || !inspection.limits.is_empty(),
        inspection,
        findings: fresh.findings,
        findings_total: fresh.findings_total,
        findings_truncated: fresh.findings_truncated,
    };
    report.validate_binding(&report.context_revision)?;
    Ok(report)
}
impl DigestPreview {
    pub fn selected_fragment(
        &self,
        choice: &crate::digest_cut_selection::Choice,
    ) -> Result<&crate::restriction_inspection::FragmentSummary, Error> {
        self.validate_binding(&self.context_revision)?;
        crate::digest_cut_selection::selected_fragment_raw_with_source(
            &self.inspection,
            choice,
            &self.source_duplex,
        )
    }
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if self.engine_version != ENGINE_VERSION
            || !hex(expected)
            || self.context_revision != expected
            || self.findings.len() != self.findings_total.min(128)
            || self.findings.iter().any(|s| s.chars().count() > 2048)
            || (self.findings_total > self.findings.len() && !self.findings_truncated)
            || self.findings_total > crate::run_lineage::MAX_RUN_BYTES
        {
            return Err(bad("digest preview binding/bounds differ"));
        }
        self.selection.validate()?;
        self.inspection.validate_with_source(&self.source_duplex)?;
        if self.inspection.enzymes
            != self
                .selection
                .enzymes
                .iter()
                .map(|pin| pin.enzyme.clone())
                .collect::<Vec<_>>()
        {
            return Err(bad(
                "Inspection enzyme geometry differs from its pinned catalog selection.",
            ));
        }
        summary_shape(&self.source_modifications)?;
        crate::digest_cut_selection::source_geometry_preflight(
            &self.inspection,
            &self.source_duplex,
        )?;
        if self.source_modifications.top_length != self.inspection.source_length
            || self.source_modifications.bottom_length != self.source_duplex.bottom_length
        {
            return Err(bad(
                "Source chemistry summaries differ from the actual source geometry.",
            ));
        }
        if self.fragment_modifications.len()
            != if self.inspection.products_determined {
                self.inspection.fragments.len()
            } else {
                0
            }
        {
            return Err(bad(
                "fragment chemistry summaries differ from actual inspection",
            ));
        }
        for (s, f) in self
            .fragment_modifications
            .iter()
            .zip(&self.inspection.fragments)
        {
            summary_shape(s)?;
            if s.top_length != f.duplex.top_length
                || s.bottom_length != f.duplex.bottom_length
                || s.generation != self.source_modifications.generation
            {
                return Err(bad("fragment summary coordinate lengths differ"));
            }
        }
        encoded(self, MAX_RECORD_BYTES)?;
        Ok(())
    }
}

// Native source-bound previews share the same origin replay and mutation kernels
// as installation. They are documentary objects, never deserialized material.
dto! {pub struct MismatchPreview{pub engine_version:String,pub context_revision:String,pub source_revision:String,
#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]pub preview:verification::MutagenesisPreview,
pub source_duplex:chemistry::DuplexSummary,pub source_modifications:ModificationSummary,
pub incomplete:bool,pub findings:Vec<String>,pub findings_total:usize,pub findings_truncated:bool}}
dto! {pub struct RepairPreview{pub engine_version:String,pub context_revision:String,pub source_revision:String,pub molecule:Molecule,
#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]pub repair:chemistry::EndRepairReport,
pub source_duplex:chemistry::DuplexSummary,pub source_modifications:ModificationSummary,
pub incomplete:bool,pub findings:Vec<String>,pub findings_total:usize,pub findings_truncated:bool}}
fn source_origin_preflight(source: &Molecule, origin: &SourceOrigin) -> Result<(), Error> {
    crate::model::bounds::model_preflight(&[source])?;
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
pub fn mismatch_pcr_preview_revision(
    source: &Molecule,
    origin: &SourceOrigin,
    forward: &Primer,
    reverse: &Primer,
    options: &crate::mutagenesis::PcrOptions,
) -> Result<String, Error> {
    crate::assembly_recipe::mismatch_pcr_shape(forward, reverse, options)?;
    source_origin_preflight(source, origin)?;
    if crate::assembly_recipe::mismatch_pcr_work(source.sequence.len(), forward, reverse, options)?
        > 200_000_000
    {
        return Err(Error::TooLarge);
    }
    crate::assembly_recipe::hash(&(
        ENGINE_VERSION,
        "mismatch-pcr-preview",
        verification::ENGINE_VERSION,
        source,
        origin,
        forward,
        reverse,
        options,
    ))
}
fn mismatch_from_fresh(
    source: &Molecule,
    origin: &SourceOrigin,
    fresh: &FreshMaterial,
    forward: &Primer,
    reverse: &Primer,
    options: &crate::mutagenesis::PcrOptions,
) -> Result<MismatchPreview, Error> {
    let context_revision =
        mismatch_pcr_preview_revision(source, origin, forward, reverse, options)?;
    let report = crate::mutagenesis::simulate_pcr_reconstructed(
        source,
        &fresh.tracked.duplex,
        forward,
        reverse,
        options,
    )?;
    let mut preview =
        verification::project_mutagenesis_preview(source, forward, reverse, options, report)?;
    for document in &mut preview.documents {
        document.metadata.insert("duplex_model".into(),"Intentional-mismatch PCR reference preview only; replay the exact native derivation to establish reusable material.".into());
        document.metadata.insert(
            "source:material_derivation:preview_context".into(),
            context_revision.clone(),
        );
    }
    let local = preview
        .findings
        .iter()
        .chain(&preview.report.warnings)
        .cloned()
        .collect::<Vec<_>>();
    let total = fresh
        .findings_total
        .checked_add(local.len())
        .ok_or(Error::TooLarge)?;
    let (findings, findings_total, findings_truncated) = limitations(
        local
            .into_iter()
            .chain(fresh.findings.iter().cloned())
            .collect(),
        total,
        fresh.findings_truncated,
    );
    let result = MismatchPreview {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        source_revision: source.revision(),
        source_duplex: fresh.tracked.duplex.summary()?,
        source_modifications: modified_duplex::summary(&fresh.tracked)?,
        incomplete: fresh.incomplete || !preview.report.enumeration_complete,
        preview,
        findings,
        findings_total,
        findings_truncated,
    };
    result.validate_binding(&result.context_revision)?;
    Ok(result)
}
pub fn preview_mismatch_pcr(
    source: &Molecule,
    origin: &SourceOrigin,
    forward: &Primer,
    reverse: &Primer,
    options: &crate::mutagenesis::PcrOptions,
) -> Result<MismatchPreview, Error> {
    mismatch_pcr_preview_revision(source, origin, forward, reverse, options)?;
    let fresh = fresh_source(source, origin)?;
    mismatch_from_fresh(source, origin, &fresh, forward, reverse, options)
}
pub fn repair_preview_revision(
    source: &Molecule,
    origin: &SourceOrigin,
    method: chemistry::EndRepairMethod,
) -> Result<String, Error> {
    source_origin_preflight(source, origin)?;
    crate::assembly_recipe::hash(&(
        ENGINE_VERSION,
        "end-repair-preview",
        verification::ENGINE_VERSION,
        source,
        origin,
        method,
    ))
}
pub fn preview_repair(
    source: &Molecule,
    origin: &SourceOrigin,
    method: chemistry::EndRepairMethod,
) -> Result<RepairPreview, Error> {
    let context_revision = repair_preview_revision(source, origin, method)?;
    let fresh = fresh_source(source, origin)?;
    let repair = chemistry::repair_ends(&fresh.tracked.duplex, method)?;
    let mut local = repair.warnings.clone();
    let mut molecule = verification::repaired_document(
        source,
        &repair,
        format!("{}_end_repair", source.name),
        &mut local,
    )?;
    molecule.metadata.insert(
        "source:material_derivation:preview_context".into(),
        context_revision.clone(),
    );
    local.push("Reference preview of explicit template-supported end repair. It restores no material authority, infers no phosphorylation or nick sealing, and checks no intended target.".into());
    let total = fresh
        .findings_total
        .checked_add(local.len())
        .ok_or(Error::TooLarge)?;
    let (findings, findings_total, findings_truncated) = limitations(
        local
            .into_iter()
            .chain(fresh.findings.iter().cloned())
            .collect(),
        total,
        fresh.findings_truncated,
    );
    let result = RepairPreview {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        source_revision: source.revision(),
        source_duplex: fresh.tracked.duplex.summary()?,
        source_modifications: modified_duplex::summary(&fresh.tracked)?,
        molecule,
        repair,
        incomplete: fresh.incomplete,
        findings,
        findings_total,
        findings_truncated,
    };
    result.validate_binding(&result.context_revision)?;
    Ok(result)
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
        || total > crate::run_lineage::MAX_RUN_BYTES
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
        || duplex.bottom_start.unsigned_abs() > crate::model::MAX_SEQUENCE_LEN as u64
        || duplex.ends.as_ref().is_some_and(|(a, b)| {
            a.overhang.len() > crate::model::MAX_SEQUENCE_LEN
                || b.overhang.len() > crate::model::MAX_SEQUENCE_LEN
        })
    {
        return Err(bad(
            "processing preview source geometry/modifications differ",
        ));
    }
    Ok(())
}
impl MismatchPreview {
    /// Exact same eligibility as installation. Incomplete enumeration stays conditional;
    /// two reported geometries never establish a selected isolated product.
    pub fn reusable(&self) -> bool {
        self.validate_binding(&self.context_revision).is_ok()
            && self.preview.report.products.len() == 1
            && self.preview.documents.len() == 1
            && self.preview.report.products[0].certain
            && !self.preview.report.products[0].overlapping_domains
    }
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        preview_receipt_shape(
            &self.engine_version,
            &self.context_revision,
            expected,
            &self.source_revision,
            &self.findings,
            self.findings_total,
            self.findings_truncated,
        )?;
        preview_source_shape(&self.source_duplex, &self.source_modifications)?;
        let report = &self.preview.report;
        if self.source_duplex.bottom_start != 0
            || self.source_duplex.top_length != self.source_duplex.bottom_length
            || !self.source_duplex.nicks.is_empty()
        {
            return Err(bad(
                "mismatch preview requires actual fully paired intact source geometry",
            ));
        }
        for binding in &report.bindings {
            if binding.source_revision != self.source_revision
                || crate::assembly_recipe::hash(&binding.options)?
                    != crate::assembly_recipe::hash(&report.options.binding)?
            {
                return Err(bad(
                    "mismatch preview binding report differs from declared policy/source",
                ));
            }
        }
        if report.source_revision != self.source_revision
            || report.products.len() != self.preview.documents.len()
            || report.products.len() > 10000
            || report.pair_checks > report.options.max_pair_checks
            || (!report.enumeration_complete && !self.incomplete)
        {
            return Err(bad("mismatch preview report correspondence differs"));
        }
        for (document, geometry) in self.preview.documents.iter().zip(&report.products) {
            crate::model::bounds::model_preflight(&[document])?;
            document.validate()?;
            if document.sequence != geometry.sequence
                || document.topology != crate::Topology::Linear
                || document
                    .metadata
                    .get("source:material_derivation:preview_context")
                    .map(String::as_str)
                    != Some(expected)
            {
                return Err(bad(
                    "mismatch preview geometry/document correspondence differs",
                ));
            }
        }
        encoded(self, 32 * 1024 * 1024)?;
        Ok(())
    }
}
impl RepairPreview {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        preview_receipt_shape(
            &self.engine_version,
            &self.context_revision,
            expected,
            &self.source_revision,
            &self.findings,
            self.findings_total,
            self.findings_truncated,
        )?;
        preview_source_shape(&self.source_duplex, &self.source_modifications)?;
        self.repair.source.validate()?;
        self.repair.duplex.validate()?;
        self.molecule.validate()?;
        if self.repair.source.summary()?.top_length != self.source_duplex.top_length
            || crate::assembly_recipe::hash(&self.repair.source.summary()?)?
                != crate::assembly_recipe::hash(&self.source_duplex)?
            || self.molecule.sequence != self.repair.duplex.top
            || self.molecule.topology != self.repair.duplex.topology
            || self
                .molecule
                .metadata
                .get("source:material_derivation:preview_context")
                .map(String::as_str)
                != Some(expected)
        {
            return Err(bad(
                "repair preview native geometry/reference correspondence differs",
            ));
        }
        encoded(self, 32 * 1024 * 1024)?;
        Ok(())
    }
}

// A bounded source-backed OE reference preview. Its documents and geometry are
// documentary; only native recipe installation/replay establishes material.
dto! {pub struct OverlapExtensionPreview{
pub engine_version:String,pub context_revision:String,pub source_revision:String,
pub source:EventBinding,pub origin_revision:Option<String>,pub origin_incomplete:bool,
pub preparation:crate::overlap_extension_pcr::Request,
pub pcr_stages:Vec<crate::overlap_extension_digest_inspection::PcrStage>,
#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]
pub extension:Option<crate::overlap_extension::Geometry>,pub product:Option<Molecule>,
pub incomplete:bool,pub findings:Vec<String>,pub findings_total:usize,pub findings_truncated:bool}}
fn overlap_extension_context(
    source: &EventBinding,
    origin: &SourceOrigin,
    preparation: &crate::overlap_extension_pcr::Request,
) -> Result<String, Error> {
    overlap_extension_context_fields(
        source,
        &match origin {
            SourceOrigin::OriginalReference => None,
            SourceOrigin::Saved(origin) => Some(crate::assembly_recipe::hash(origin)?),
        },
        preparation,
    )
}
fn overlap_extension_context_fields(
    source: &EventBinding,
    origin_revision: &Option<String>,
    preparation: &crate::overlap_extension_pcr::Request,
) -> Result<String, Error> {
    crate::assembly_recipe::hash(&(
        ENGINE_VERSION,
        "overlap-extension-preview",
        verification::ENGINE_VERSION,
        source,
        origin_revision,
        preparation,
    ))
}
/// Admission is workspace-bound even for original references: content equality
/// cannot substitute for the source event selected by the user.
pub fn overlap_extension_preview_revision(
    w: &Workspace,
    source: &EventBinding,
    origin: &SourceOrigin,
    preparation: &crate::overlap_extension_pcr::Request,
) -> Result<String, Error> {
    preparation.validate_shape()?;
    let molecule = workspace_source(w, source, origin)?;
    if matches!(origin, SourceOrigin::OriginalReference)
        && [&preparation.left, &preparation.right]
            .iter()
            .any(|pair| pair.template_strand != chemistry::Strand::Top)
    {
        return Err(bad(
            "Original-reference OE PCR requires Top source strands; Bottom requires native source replay",
        ));
    }
    crate::overlap_extension_execution::preflight(molecule, preparation)?;
    overlap_extension_context(source, origin, preparation)
}
pub fn preview_overlap_extension(
    w: &Workspace,
    source: &EventBinding,
    origin: &SourceOrigin,
    preparation: &crate::overlap_extension_pcr::Request,
) -> Result<OverlapExtensionPreview, Error> {
    let context_revision = overlap_extension_preview_revision(w, source, origin, preparation)?;
    let molecule = run_material::resolve(w, source)?;
    let fresh = fresh_source(molecule, origin)?;
    let executed =
        crate::overlap_extension_execution::execute(molecule, Some(&fresh.tracked), preparation)?;
    let total = fresh
        .findings_total
        .checked_add(executed.findings.len())
        .ok_or(Error::TooLarge)?;
    let (findings, findings_total, findings_truncated) = limitations(
        executed
            .findings
            .into_iter()
            .chain(fresh.findings)
            .collect(),
        total,
        fresh.findings_truncated,
    );
    let result = OverlapExtensionPreview {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        source_revision: molecule.revision(),
        source: source.clone(),
        origin_revision: match origin {
            SourceOrigin::OriginalReference => None,
            SourceOrigin::Saved(origin) => Some(crate::assembly_recipe::hash(origin)?),
        },
        origin_incomplete: fresh.incomplete,
        preparation: preparation.clone(),
        pcr_stages: executed.stages,
        extension: executed.extension,
        product: executed.product,
        incomplete: fresh.incomplete || executed.incomplete,
        findings,
        findings_total,
        findings_truncated,
    };
    result.validate_binding(&result.context_revision)?;
    Ok(result)
}
impl OverlapExtensionPreview {
    /// No source replay, PCR search or target evaluation occurs here.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.preparation.validate_shape()?;
        run_material::event_shape(&self.source)?;
        if self.engine_version != ENGINE_VERSION
            || !hex(expected)
            || self.context_revision != expected
            || self.source_revision != self.source.revision
            || self.origin_revision.as_ref().is_some_and(|s| !hex(s))
            || self.origin_revision.is_none() && self.origin_incomplete
            || self.context_revision
                != overlap_extension_context_fields(
                    &self.source,
                    &self.origin_revision,
                    &self.preparation,
                )?
            || self.findings.len() != self.findings_total.min(128)
            || self.findings_total > crate::run_lineage::MAX_RUN_BYTES
            || self.findings.iter().any(|s| s.chars().count() > 2048)
            || self.findings_total > self.findings.len() && !self.findings_truncated
        {
            return Err(bad(
                "OE preview source, origin, preparation or bounded context changed",
            ));
        }
        crate::overlap_extension_execution::validate_evidence(
            &self.source_revision,
            &self.preparation,
            &self.pcr_stages,
            &self.extension,
            &self.product,
        )?;
        let expected_incomplete = self.origin_incomplete
            || self.product.is_none()
            || self
                .extension
                .as_ref()
                .is_none_or(|g| !g.terminal_pairings.complete || !g.terminal_pairings.unique);
        if self.incomplete != expected_incomplete {
            return Err(bad(
                "OE preview lost inherited or terminal-pair uncertainty",
            ));
        }
        encoded(self, crate::run_lineage::MAX_RUN_BYTES)?;
        Ok(())
    }
    /// Conditional source/terminal uncertainty is retained, never upgraded by selection.
    pub fn reusable(&self) -> bool {
        self.product.is_some()
            && self.extension.is_some()
            && self.pcr_stages.len() == 3
            && self
                .pcr_stages
                .iter()
                .all(|p| p.geometries.len() == 1 && p.product_revision.is_some())
    }
    pub fn selection(&self) -> Result<OverlapExtensionSelection, Error> {
        self.validate_binding(&self.context_revision)?;
        if !self.reusable() {
            return Err(bad(
                "No reusable OE outer-PCR product: each source PCR and outer PCR must have one actual geometry",
            ));
        }
        Ok(OverlapExtensionSelection {
            expected_preview_revision: self.context_revision.clone(),
            expected_product_revision: self
                .product
                .as_ref()
                .ok_or_else(|| bad("outer-PCR product absent"))?
                .revision(),
        })
    }
}

// Documentary source/chemistry comparison; serialized summaries never restore state.
dto! {pub struct TerminalProcessingPreview{
pub engine_version:String,pub context_revision:String,pub source_revision:String,pub source_sequence_blake3:String,
pub source:EventBinding,pub origin_revision:Option<String>,pub origin_incomplete:bool,pub action:TerminalAction,
pub product:Molecule,pub source_duplex:chemistry::DuplexSummary,pub duplex:chemistry::DuplexSummary,
pub source_modifications:ModificationSummary,pub modifications:ModificationSummary,
pub incomplete:bool,pub findings:Vec<String>,pub findings_total:usize,pub findings_truncated:bool}}
fn source_sequence_hash(source: &Molecule) -> String {
    blake3::hash(source.sequence.as_str().as_bytes())
        .to_hex()
        .to_string()
}
fn workspace_source<'a>(
    w: &'a Workspace,
    source: &EventBinding,
    origin: &SourceOrigin,
) -> Result<&'a Molecule, Error> {
    origin.validate_shallow_for(source)?;
    crate::run_material::material_origins_preflight(match origin {
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
fn terminal_context_fields(
    source: &EventBinding,
    source_sequence_blake3: &str,
    origin_revision: &Option<String>,
    action: &TerminalAction,
) -> Result<String, Error> {
    crate::assembly_recipe::hash(&(
        ENGINE_VERSION,
        "terminal-processing-preview",
        verification::ENGINE_VERSION,
        source,
        source_sequence_blake3,
        origin_revision,
        action,
    ))
}
pub fn terminal_processing_preview_revision(
    w: &Workspace,
    source: &EventBinding,
    origin: &SourceOrigin,
    action: &TerminalAction,
) -> Result<String, Error> {
    let molecule = workspace_source(w, source, origin)?;
    terminal_context_fields(
        source,
        &source_sequence_hash(molecule),
        &match origin {
            SourceOrigin::OriginalReference => None,
            SourceOrigin::Saved(origin) => Some(crate::assembly_recipe::hash(origin)?),
        },
        action,
    )
}
fn terminal_from_fresh(
    source: &Molecule,
    fresh: &FreshMaterial,
    action: &TerminalAction,
    name: &str,
) -> Result<(Molecule, TrackedDuplex, Vec<String>), Error> {
    let (tracked, findings) = match action {
        TerminalAction::TreatEnds { treatment } => {
            let tracked = fresh.tracked.treat(*treatment)?;
            (
                tracked,
                vec![format!(
                    "Explicit {treatment:?} changes only the two external 5′ groups on actual intact linear strands. Strand letters, 3′ groups and positional modifications remain; no enzyme yield, nick repair or independent-target success is inferred."
                )],
            )
        }
        TerminalAction::Tail3Prime { base } => {
            let tracked = fresh.tracked.tail(*base)?;
            (tracked,vec![format!("Exactly one explicitly declared {} base is added to each actual 3′ hydroxyl on the blunt linear duplex. Source calls, 5′ groups, existing nicks and modifications remain; new calls are unmodified under the existing standard-nucleotide model.",base.as_str()),"This ideal endpoint transform does not predict enzyme/kit activity, tail heterogeneity, yield or experimental completion. Symmetric single-base tails do not establish insertion orientation or unique products; later ligation must retain alternative products and the independent intended path.".into()])
        }
    };
    let mut reference = tracked.duplex.molecule(name)?;
    // Source top-reference calls remain in precisely the same coordinates. Appended
    // terminal bases are outside every original feature, never inserted into a CDS.
    reference.features = source.features.clone();
    reference.metadata.insert(
        "source:material_derivation:terminal_source_revision".into(),
        source.revision(),
    );
    reference.metadata.insert(
        "source:material_derivation:terminal_source_sequence_blake3".into(),
        source_sequence_hash(source),
    );
    reference.metadata.insert(
        "source:material_derivation:terminal_action".into(),
        serde_json::to_string(action).map_err(|_| Error::TooLarge)?,
    );
    reference.metadata.insert("source:material_derivation:terminal_scope".into(),"Native source-backed endpoint processing. Detached reference and summaries are not reusable stock; retain exact source events and replay the native recipe. No intended target was checked.".into());
    reference.metadata.insert("source:material_derivation:scope".into(),"Terminal processing reference preview. Exact native source ancestry and processing recipe are required for material reuse.".into());
    reference.validate()?;
    Ok((reference, tracked, findings))
}
pub fn preview_terminal_processing(
    w: &Workspace,
    source: &EventBinding,
    origin: &SourceOrigin,
    action: &TerminalAction,
) -> Result<TerminalProcessingPreview, Error> {
    let context_revision = terminal_processing_preview_revision(w, source, origin, action)?;
    let molecule = run_material::resolve(w, source)?;
    let fresh = fresh_source(molecule, origin)?;
    let (product, tracked, local) =
        terminal_from_fresh(molecule, &fresh, action, "Terminally processed material")?;
    let total = fresh
        .findings_total
        .checked_add(local.len())
        .ok_or(Error::TooLarge)?;
    let (findings, findings_total, findings_truncated) = limitations(
        local
            .into_iter()
            .chain(fresh.findings.iter().cloned())
            .collect(),
        total,
        fresh.findings_truncated,
    );
    let result = TerminalProcessingPreview {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        source_revision: molecule.revision(),
        source_sequence_blake3: source_sequence_hash(molecule),
        source: source.clone(),
        origin_revision: match origin {
            SourceOrigin::OriginalReference => None,
            SourceOrigin::Saved(origin) => Some(crate::assembly_recipe::hash(origin)?),
        },
        origin_incomplete: fresh.incomplete,
        action: *action,
        product,
        source_duplex: fresh.tracked.duplex.summary()?,
        duplex: tracked.duplex.summary()?,
        source_modifications: modified_duplex::summary(&fresh.tracked)?,
        modifications: modified_duplex::summary(&tracked)?,
        incomplete: fresh.incomplete,
        findings,
        findings_total,
        findings_truncated,
    };
    result.validate_binding(&result.context_revision)?;
    Ok(result)
}
impl TerminalProcessingPreview {
    /// Bounded documentary correspondence only. Never replays source material or
    /// executes endpoint processing, and cannot authorize installation on its own.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        run_material::event_shape(&self.source)?;
        crate::model::bounds::model_preflight(&[&self.product])?;
        self.product.validate()?;
        summary_shape(&self.source_modifications)?;
        summary_shape(&self.modifications)?;
        if self.engine_version != ENGINE_VERSION
            || !hex(expected)
            || self.context_revision != expected
            || self.source_revision != self.source.revision
            || !hex(&self.source_sequence_blake3)
            || self.origin_revision.as_ref().is_some_and(|s| !hex(s))
            || self.origin_revision.is_none() && self.origin_incomplete
            || self.context_revision
                != terminal_context_fields(
                    &self.source,
                    &self.source_sequence_blake3,
                    &self.origin_revision,
                    &self.action,
                )?
            || self.incomplete != self.origin_incomplete
            || self.findings.len() != self.findings_total.min(128)
            || self.findings_total > crate::run_lineage::MAX_RUN_BYTES
            || self.findings.iter().any(|s| s.chars().count() > 2048)
            || self.findings_total > self.findings.len() && !self.findings_truncated
            || self.source_duplex.top_length != self.source_modifications.top_length
            || self.source_duplex.bottom_length != self.source_modifications.bottom_length
            || self.duplex.top_length != self.modifications.top_length
            || self.duplex.bottom_length != self.modifications.bottom_length
            || self.product.sequence.alphabet() != Alphabet::Dna
            || self.product.sequence.len() != self.duplex.top_length
            || self.product.topology != self.duplex.topology
            || self.duplex.topology != crate::Topology::Linear
            || self.source_duplex.topology != crate::Topology::Linear
            || self.source_duplex.nicks != self.duplex.nicks
            || self.product.metadata.get("duplex_model")
                != Some(&serde_json::to_string(&self.duplex).map_err(|_| Error::TooLarge)?)
            || self
                .product
                .metadata
                .get("source:material_derivation:terminal_source_revision")
                != Some(&self.source_revision)
            || self
                .product
                .metadata
                .get("source:material_derivation:terminal_source_sequence_blake3")
                != Some(&self.source_sequence_blake3)
            || self
                .product
                .metadata
                .get("source:material_derivation:terminal_action")
                != Some(&serde_json::to_string(&self.action).map_err(|_| Error::TooLarge)?)
            || self.source_modifications.known_modified != self.modifications.known_modified
            || self.source_modifications.unknown_eligible != self.modifications.unknown_eligible
            || self.source_modifications.generation != self.modifications.generation
        {
            return Err(bad(
                "Terminal preview changed its exact source, action, chemistry, inherited limitations or bounded findings",
            ));
        }
        let source_ends = self
            .source_duplex
            .ends
            .as_ref()
            .ok_or_else(|| bad("Terminal processing source has no free ends"))?;
        let product_ends = self
            .duplex
            .ends
            .as_ref()
            .ok_or_else(|| bad("Terminal processing product has no free ends"))?;
        match self.action {
            TerminalAction::TreatEnds { treatment } => {
                if !self.duplex.nicks.is_empty()
                    || self.duplex.top_length != self.source_duplex.top_length
                    || self.duplex.bottom_length != self.source_duplex.bottom_length
                    || self.duplex.bottom_start != self.source_duplex.bottom_start
                    || source_sequence_hash(&self.product) != self.source_sequence_blake3
                    || self.modifications.known_unmodified
                        != self.source_modifications.known_unmodified
                {
                    return Err(bad(
                        "Terminal treatment changed source geometry, letters or modifications",
                    ));
                }
                for (before, after) in [
                    (&source_ends.0, &product_ends.0),
                    (&source_ends.1, &product_ends.1),
                ] {
                    let mut wanted = before.clone();
                    wanted.chemistry.five_prime = match treatment {
                        chemistry::EndTreatment::Phosphorylate => chemistry::FivePrime::Phosphate,
                        chemistry::EndTreatment::Dephosphorylate => chemistry::FivePrime::Hydroxyl,
                    };
                    if &wanted != after {
                        return Err(bad(
                            "Terminal treatment changed more than the declared external 5′ groups",
                        ));
                    }
                }
            }
            TerminalAction::Tail3Prime { base } => {
                let length = self.source_duplex.top_length;
                if self.source_duplex.bottom_start != 0
                    || self.source_duplex.bottom_length != length
                    || self.duplex.bottom_start != -1
                    || self.duplex.top_length != length.checked_add(1).ok_or(Error::TooLarge)?
                    || self.duplex.bottom_length != self.duplex.top_length
                    || self.product.sequence.as_str().as_bytes().last().copied()
                        != base.as_str().as_bytes().first().copied()
                    || blake3::hash(&self.product.sequence.as_str().as_bytes()[..length])
                        .to_hex()
                        .as_str()
                        != self.source_sequence_blake3
                    || self.modifications.known_unmodified
                        != self
                            .source_modifications
                            .known_unmodified
                            .checked_add(
                                if matches!(
                                    base,
                                    crate::tailing::TailingBase::A | crate::tailing::TailingBase::C
                                ) {
                                    2
                                } else {
                                    0
                                },
                            )
                            .ok_or(Error::TooLarge)?
                {
                    return Err(bad(
                        "Terminal tailing differs from exactly one declared newly synthesized base per strand",
                    ));
                }
                for (before, after) in [
                    (&source_ends.0, &product_ends.0),
                    (&source_ends.1, &product_ends.1),
                ] {
                    if before.polarity != chemistry::Polarity::Blunt
                        || !before.overhang.is_empty()
                        || before.chemistry.three_prime != chemistry::ThreePrime::Hydroxyl
                        || after.polarity != chemistry::Polarity::ThreePrime
                        || after.overhang != base.as_str()
                        || after.chemistry != before.chemistry
                    {
                        return Err(bad(
                            "Terminal tailing changed inherited groups or did not expose the declared 3′ tails",
                        ));
                    }
                }
            }
        }
        encoded(self, crate::run_lineage::MAX_RUN_BYTES)?;
        Ok(())
    }
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
    crate::model::bounds::model_preflight(&[molecule])?;
    verification::declared_oligo(molecule)?;
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
        crate::model::bounds::model_preflight(&[&snapshot.molecule])?;
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
fn anneal_options_shape(options: &chemistry::AnnealOptions) -> Result<(), Error> {
    crate::assembly_recipe::anneal_shape("bottom", options, 0)
}
fn candidate_id_shape(id: &str) -> bool {
    id.strip_prefix("anneal-").is_some_and(hex)
}
fn anneal_context(
    sources: &[SourceBinding],
    options: &chemistry::AnnealOptions,
) -> Result<String, Error> {
    crate::assembly_recipe::hash(&(
        ENGINE_VERSION,
        "annealed-oligo-preview",
        verification::ENGINE_VERSION,
        sources,
        options,
    ))
}
fn fresh(
    snapshots: &[SourceSnapshot],
    request: &Request,
    revision: &str,
) -> Result<FreshMaterial, Error> {
    validate_source_snapshots(request, snapshots)?;
    if let Derivation::AnnealedOligos { options, selection } = &request.derivation {
        let mut findings = Vec::new();
        let preview = verification::simulate_oligo_annealing(
            &snapshots[0].molecule,
            &snapshots[1].molecule,
            options,
        )?;
        if anneal_context(&request.source_bindings(), options)?
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
        let mut reference = verification::annealed_document_inner(
            &snapshots[0].molecule,
            &snapshots[1].molecule,
            &preview.annealing,
            index,
            &mut findings,
        )?;
        protect_anneal_preview(&mut reference)?;
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
            reference,
            tracked,
            incomplete,
            findings,
            findings_total,
            findings_truncated,
        })
    } else {
        fresh_template(&snapshots[0].molecule, request, revision)
    }
}
fn protect_anneal_preview(product: &mut Molecule) -> Result<(), Error> {
    product.metadata.insert("source:material_derivation:scope".into(),"Source-backed standard-oligo annealing reference preview. Both exact native source events and the selected geometry recipe are required for reusable material.".into());
    product.validate()
}
dto! {pub struct AnnealedOligosPreview{
pub engine_version:String,pub context_revision:String,pub sources:Vec<SourceBinding>,pub options:chemistry::AnnealOptions,
#[cfg_attr(feature="json-schema",schemars(with="serde_json::Value"))]pub preview:verification::AnnealPreview,
pub products:Vec<Molecule>,pub incomplete:bool,pub findings:Vec<String>}}
pub fn annealed_oligos_preview_revision(
    w: &Workspace,
    top: &EventBinding,
    bottom: &EventBinding,
    options: &chemistry::AnnealOptions,
) -> Result<String, Error> {
    anneal_options_shape(options)?;
    run_material::event_shape(top)?;
    run_material::event_shape(bottom)?;
    w.validate()?;
    for event in [top, bottom] {
        declared_oligo_source(run_material::resolve(w, event)?)?;
    }
    anneal_context(
        &[
            SourceBinding {
                role: SourceRole::TopOligo,
                event: top.clone(),
            },
            SourceBinding {
                role: SourceRole::BottomOligo,
                event: bottom.clone(),
            },
        ],
        options,
    )
}
pub fn preview_annealed_oligos(
    w: &Workspace,
    top: &EventBinding,
    bottom: &EventBinding,
    options: &chemistry::AnnealOptions,
) -> Result<AnnealedOligosPreview, Error> {
    let context_revision = annealed_oligos_preview_revision(w, top, bottom, options)?;
    let a = run_material::resolve(w, top)?;
    let b = run_material::resolve(w, bottom)?;
    let preview = verification::simulate_oligo_annealing(a, b, options)?;
    let mut findings = preview.annealing.warnings.clone();
    let mut products = Vec::new();
    for index in 0..preview.annealing.candidates.len() {
        let mut product =
            verification::annealed_document_inner(a, b, &preview.annealing, index, &mut findings)?;
        protect_anneal_preview(&mut product)?;
        products.push(product);
    }
    let result = AnnealedOligosPreview {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        sources: vec![
            SourceBinding {
                role: SourceRole::TopOligo,
                event: top.clone(),
            },
            SourceBinding {
                role: SourceRole::BottomOligo,
                event: bottom.clone(),
            },
        ],
        options: options.clone(),
        incomplete: !preview.annealing.complete || !preview.annealing.unique,
        preview,
        products,
        findings,
    };
    result.validate_binding(&result.context_revision)?;
    Ok(result)
}
impl AnnealedOligosPreview {
    /// Documentary binding only: does not rerun annealing or establish private state.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        validate_source_bindings(&self.sources)?;
        anneal_options_shape(&self.options)?;
        if self.sources.len() != 2
            || self.engine_version != ENGINE_VERSION
            || !hex(expected)
            || self.context_revision != expected
            || self.context_revision != anneal_context(&self.sources, &self.options)?
            || self.preview.top_revision != self.sources[0].event.revision
            || self.preview.bottom_revision != self.sources[1].event.revision
            || self.preview.annealing.options != self.options
            || self.preview.annealing.engine_version != chemistry::ANNEAL_ENGINE_VERSION
            || self.products.len() != self.preview.annealing.candidates.len()
            || self.products.len() > self.options.max_candidates
            || self.preview.annealing.comparisons > self.options.max_comparisons
            || self.preview.annealing.unique
                != (self.preview.annealing.complete && self.products.len() == 1)
            || self.incomplete
                != (!self.preview.annealing.complete || !self.preview.annealing.unique)
            || self.findings.len() > 1000
            || self.findings.iter().any(|s| s.len() > 8192)
        {
            return Err(bad(
                "Annealed preview changed source roles, candidate evidence or context",
            ));
        }
        let mut ids = BTreeSet::new();
        for (candidate, product) in self.preview.annealing.candidates.iter().zip(&self.products) {
            candidate.duplex.validate()?;
            product.validate()?;
            if !candidate_id_shape(&candidate.id)
                || !ids.insert(&candidate.id)
                || product.sequence != candidate.duplex.top
                || product.topology != crate::Topology::Linear
                || product.metadata.get("anneal:top_revision")
                    != Some(&self.sources[0].event.revision)
                || product.metadata.get("anneal:bottom_revision")
                    != Some(&self.sources[1].event.revision)
                || product.metadata.get("anneal:candidate_id") != Some(&candidate.id)
                || product.metadata.get("duplex_model")
                    != Some(
                        &serde_json::to_string(&candidate.duplex.summary()?)
                            .map_err(|_| Error::TooLarge)?,
                    )
                || !product
                    .metadata
                    .contains_key("source:material_derivation:scope")
            {
                return Err(bad(
                    "Annealed product differs from its retained candidate/source geometry",
                ));
            }
        }
        encoded(self, crate::run_lineage::MAX_RUN_BYTES)?;
        Ok(())
    }
    pub fn selection(&self, candidate_id: &str) -> Result<AnnealSelection, Error> {
        self.validate_binding(&self.context_revision)?;
        let index = self
            .preview
            .annealing
            .candidates
            .iter()
            .position(|c| c.id == candidate_id)
            .ok_or_else(|| bad("choose a retained native annealing candidate ID"))?;
        Ok(AnnealSelection {
            expected_preview_revision: self.context_revision.clone(),
            candidate_id: candidate_id.into(),
            expected_product_revision: self.products[index].revision(),
        })
    }
}
