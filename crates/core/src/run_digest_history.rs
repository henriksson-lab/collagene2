//! Target-free saved-run digest derivations. Records and capsules contain recipes,
//! never serialized physical state or an invented independent target assessment.
use crate::{
    Error, Molecule, assembly_recipe,
    catalog_digest::{self, CatalogSelection},
    chemistry,
    modified_duplex::{self, TrackedDuplex},
    restriction_history,
    restriction_inspection::{self, InspectionOptions},
    run_lineage::{self, EventBinding},
    run_material, transfer, verification,
    workspace::{self, Workspace, WorkspaceAttachment},
};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};
pub const ENGINE_VERSION: &str = "saved-run-digest-history-1";
pub const WORKSPACE_ATTACHMENT_KIND: &str = "saved-run-digest-fragment-v1";
pub const MAX_RECORD_BYTES: usize = 16 << 20;
pub const MAX_CAPSULE_BYTES: usize = 8 << 20;
pub const MAX_RECORDS: usize = 128;
const CONTEXT_KEY: &str = "source:run_digest_history:context";
macro_rules! dto{($($t:tt)*)=>{#[derive(Debug,Clone,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Request{pub source:EventBinding,pub origin:run_material::Capsule,pub selection:CatalogSelection,pub options:InspectionOptions,pub fragment_index:usize}}
dto! {pub struct Record{pub schema_version:u32,pub engine_version:String,pub scientific_engine:String,pub source_workspace_fingerprint:String,pub context_revision:String,pub request:Request,pub product:EventBinding,pub product_state_revision:String,pub origin_incomplete:bool}}
dto! {pub struct Capsule{pub schema_version:u32,pub record:Record,pub source:Molecule}}
dto! {pub struct Receipt{pub engine_version:String,pub context_revision:String,pub record_hash:String,pub source:EventBinding,pub product:EventBinding,pub origin_incomplete:bool,pub origin_findings:Vec<String>,pub origin_findings_total:usize,pub origin_findings_truncated:bool,pub duplex:chemistry::DuplexSummary,pub modifications:modified_duplex::ModificationSummary}}
pub struct PreparedInstallation {
    pub workspace: Workspace,
    pub receipt: Receipt,
}
dto! {pub struct Installed{pub archive:Vec<u8>,pub receipt:Receipt}}
/// Only a successful in-process recipe replay constructs this receipt. None of
/// its fields is deserializable physical authority.
#[derive(Debug)]
pub struct FreshFragment {
    reference: Molecule,
    tracked: TrackedDuplex,
    origin_incomplete: bool,
    origin_findings: Vec<String>,
    origin_findings_total: usize,
    origin_findings_truncated: bool,
}
impl FreshFragment {
    pub fn reference(&self) -> &Molecule {
        &self.reference
    }
    pub fn tracked(&self) -> &TrackedDuplex {
        &self.tracked
    }
    pub fn origin_incomplete(&self) -> bool {
        self.origin_incomplete
    }
    pub fn origin_findings(&self) -> &[String] {
        &self.origin_findings
    }
    pub fn origin_findings_total(&self) -> usize {
        self.origin_findings_total
    }
    pub fn origin_findings_truncated(&self) -> bool {
        self.origin_findings_truncated
    }
    /// Move freshly replayed state into the verifier; never deserialize it.
    pub fn into_parts(self) -> (Molecule, TrackedDuplex, bool, Vec<String>, usize, bool) {
        (
            self.reference,
            self.tracked,
            self.origin_incomplete,
            self.origin_findings,
            self.origin_findings_total,
            self.origin_findings_truncated,
        )
    }
}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("saved-run digest history: {s}"))
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
            return Err(std::io::Error::other("digest history size bound"));
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
        bytes: vec![],
        limit,
    };
    serde_json::to_writer(&mut out, value).map_err(|_| Error::TooLarge)?;
    Ok(out.bytes)
}
fn context(fingerprint: &str, request: &Request, science: &str) -> Result<String, Error> {
    assembly_recipe::hash(&(ENGINE_VERSION, science, fingerprint, request))
}
/// Visit native plans through the common mixed-origin traversal. One nesting
/// level is reserved for this derivation; digest children are never opaque leaves.
fn origins(
    origin: &run_material::Capsule,
    mut visit: impl FnMut(&run_material::Capsule, &transfer::Candidate) -> Result<(), Error>,
) -> Result<(), Error> {
    crate::run_material_capsule::digest_origin_walk(origin, |o, c| visit(o, c), |_, _| Ok(()))
}

pub(crate) fn request_shape_shallow(request: &Request) -> Result<(), Error> {
    run_material::event_shape(&request.source)?;
    if request.origin.selection.event != request.source {
        return Err(bad("source event differs from selected native origin"));
    }
    request.selection.validate()?;
    request.options.validate()?;
    if request.fragment_index >= request.options.max_fragments {
        return Err(bad("fragment index exceeds declared inspection bounds"));
    }
    encoded(request, MAX_RECORD_BYTES)?;
    Ok(())
}
fn request_shape(request: &Request) -> Result<(), Error> {
    origins(&request.origin, |_, _| Ok(()))?;
    request_shape_shallow(request)
}
/// Exact immutable source/target/product/native-run dependencies for archive
/// closure code. These documentary pins confer no current material status.
pub fn source_bindings(record: &Record) -> Result<Vec<EventBinding>, Error> {
    record.validate_documentary()?;
    let mut found = BTreeMap::new();
    let mut add = |b: &EventBinding| {
        found.insert(
            (
                b.document_id.clone(),
                b.event_id.clone(),
                b.revision.clone(),
            ),
            b.clone(),
        );
    };
    add(&record.request.source);
    let mut nested = Vec::new();
    let mut derived = Vec::new();
    crate::run_material_capsule::digest_origin_walk_mixed(
        &record.request.origin,
        |_, _| Ok(()),
        |r, _| {
            nested.push(r.request.source.clone());
            nested.push(r.product.clone());
            Ok(())
        },
        |r, _| {
            derived.push(r.request.source.clone());
            derived.push(r.product.clone());
            Ok(())
        },
    )?;
    for binding in nested.iter().chain(&derived) {
        add(binding);
    }
    origins(&record.request.origin, |origin, _| {
        add(&origin.selection.event);
        add(&origin.target_event);
        for b in &origin.inputs {
            add(&b.event);
        }
        for p in &origin.products {
            add(&EventBinding {
                document_id: p.document_id.clone(),
                event_id: p.event_id.clone(),
                revision: p.revision.clone(),
            });
        }
        Ok(())
    })?;
    Ok(found.into_values().collect())
}
pub fn required_run_ids(record: &Record) -> Result<Vec<String>, Error> {
    record.validate_documentary()?;
    let mut ids = BTreeSet::new();
    origins(&record.request.origin, |origin, _| {
        ids.insert(origin.selection.run_id.clone());
        Ok(())
    })?;
    Ok(ids.into_iter().collect())
}
impl Record {
    pub fn validate_documentary(&self) -> Result<(), Error> {
        crate::run_material_capsule::digest_record_preflight(self)
    }
    pub(crate) fn validate_shallow(&self) -> Result<(), Error> {
        if self.schema_version != 1
            || self.engine_version != ENGINE_VERSION
            || self.scientific_engine.is_empty()
            || self.scientific_engine.len() > 128
            || !self.scientific_engine.is_ascii()
            || self.scientific_engine.bytes().any(|b| b.is_ascii_control())
            || !hex(&self.source_workspace_fingerprint)
            || !hex(&self.context_revision)
            || !hex(&self.product_state_revision)
        {
            return Err(bad("unsupported derivation schema/engine or fingerprint"));
        }
        request_shape_shallow(&self.request)?;
        run_material::event_shape(&self.product)?;
        if context(
            &self.source_workspace_fingerprint,
            &self.request,
            &self.scientific_engine,
        )? != self.context_revision
        {
            return Err(bad("derivation context differs from its exact request"));
        }
        encoded(self, MAX_RECORD_BYTES)?;
        Ok(())
    }
}
pub fn load_record(bytes: &[u8]) -> Result<Record, Error> {
    if bytes.len() > MAX_RECORD_BYTES {
        return Err(Error::TooLarge);
    }
    let record: Record = serde_json::from_slice(bytes).map_err(|e| bad(&e.to_string()))?;
    record.validate_documentary()?;
    if encoded(&record, MAX_RECORD_BYTES)? != bytes {
        return Err(bad("noncanonical derivation record"));
    }
    Ok(record)
}
fn retained_run(w: &Workspace, id: &str) -> Result<run_lineage::RunRecord, Error> {
    let mut selected = None;
    for a in w
        .attachments
        .values()
        .filter(|a| a.kind == run_lineage::WORKSPACE_ATTACHMENT_KIND)
    {
        let r = run_lineage::load_run_record(&a.bytes)?;
        if r.run_id == id {
            if selected.is_some() {
                return Err(bad("duplicate native origin run"));
            }
            selected = Some(r);
        }
    }
    selected.ok_or_else(|| bad("original native run is not retained"))
}
/// Structural retained-closure attestation. Does not call Workspace::validate or
/// scientific replay, so archive validators can use this without recursion.
fn attest(w: &Workspace, request: &Request) -> Result<(), Error> {
    request_shape_shallow(request)?;
    run_material::resolve(w, &request.source)?;
    crate::run_material_capsule::validate_digest_request_origins(w, &request.origin)
}
/// Compare one already-budgeted digest declaration to retained immutable native
/// records. Descendants are attested by the caller's SAME iterative walk.
pub(crate) fn attest_record_shallow(
    w: &Workspace,
    record: &Record,
    source: Option<&Molecule>,
) -> Result<(), Error> {
    let retained = record_for(w, &record.product)?;
    if encoded(&retained, MAX_RECORD_BYTES)? != encoded(record, MAX_RECORD_BYTES)? {
        return Err(bad(
            "digest declaration differs from retained derivation record",
        ));
    }
    let actual = run_material::resolve(w, &record.request.source)?;
    if source.is_some_and(|source| source != actual) {
        return Err(bad(
            "digest declaration source differs from retained source event",
        ));
    }
    let reference = run_material::resolve(w, &record.product)?;
    let event = w
        .history
        .iter()
        .find(|e| e.id == record.product.event_id)
        .ok_or_else(|| bad("digest event unavailable"))?;
    if event.parents != vec![record.request.source.event_id.clone()]
        || reference.metadata.get(CONTEXT_KEY) != Some(&record.context_revision)
        || !reference.metadata.contains_key("duplex_model")
    {
        return Err(bad(
            "digest product ancestry or restrictive native context changed",
        ));
    }
    Ok(())
}

fn require_current_origins(w: &Workspace, request: &Request) -> Result<(), Error> {
    origins(&request.origin, |origin, _| {
        if retained_run(w, &origin.selection.run_id)?
            .assessment
            .engine_version
            != verification::ENGINE_VERSION
        {
            return Err(bad(
                "original recorded engine is historical; review and capture a current native origin before deriving fresh material",
            ));
        }
        Ok(())
    })
}
pub fn preflight(w: &Workspace, request: &Request) -> Result<(), Error> {
    w.validate()?;
    attest(w, request)?;
    require_current_origins(w, request)?;
    let source = restriction_history::current_source(w, &request.source)?;
    let document = w
        .documents
        .iter()
        .find(|d| d.id == request.source.document_id)
        .ok_or_else(|| bad("source document missing"))?;
    if document.head != request.source.revision
        || document.head_event.as_deref() != Some(request.source.event_id.as_str())
    {
        return Err(bad("active source head/event changed"));
    }
    let c = catalog_digest::builtin_catalog()?;
    catalog_digest::source_preflight(
        &c,
        &request.selection,
        &restriction_history::source_role(&request.source)?,
        source,
        &request.options,
    )?;
    if w.attachments
        .values()
        .filter(|a| a.kind == WORKSPACE_ATTACHMENT_KIND)
        .count()
        >= MAX_RECORDS
        || w.attachments.len() >= workspace::MAX_ATTACHMENTS
    {
        return Err(Error::TooLarge);
    }
    Ok(())
}
pub fn installation_revision(w: &Workspace, request: &Request) -> Result<String, Error> {
    preflight(w, request)?;
    context(
        &run_lineage::workspace_fingerprint(w)?,
        request,
        verification::ENGINE_VERSION,
    )
}
fn fresh(
    source: &Molecule,
    request: &Request,
    context_revision: &str,
) -> Result<FreshFragment, Error> {
    request_shape(request)?;
    crate::model::bounds::model_preflight(&[source])?;
    source.validate()?;
    if source.revision() != request.source.revision {
        return Err(bad("original source model differs from exact event"));
    }
    let replay = crate::run_material_capsule::replay(&request.origin, source)?;
    let report = replay.evidence.decode_report()?;
    let origin_findings_total = report
        .findings
        .len()
        .checked_add(report.assumptions.len())
        .ok_or(Error::TooLarge)?;
    let origin_findings_truncated = origin_findings_total > 128
        || report
            .findings
            .iter()
            .chain(&report.assumptions)
            .take(128)
            .any(|s| s.chars().count() > 2048);
    let origin_findings = report
        .findings
        .iter()
        .chain(&report.assumptions)
        .take(128)
        .map(|s| s.chars().take(2048).collect())
        .collect();
    let c = catalog_digest::builtin_catalog()?;
    let enzymes = request.selection.resolved_enzymes(&c)?;
    let role = restriction_history::source_role(&request.source)?;
    let mut selected = None;
    let inspection = restriction_inspection::inspect_with_digest(
        &role,
        source,
        &enzymes,
        &request.options,
        || {
            let fragments = modified_duplex::digest(&replay.tracked, &enzymes)?;
            let mut geometries = Vec::with_capacity(fragments.len());
            for (index, fragment) in fragments.into_iter().enumerate() {
                geometries.push(chemistry::DigestFragment {
                    duplex: fragment.duplex.duplex.clone(),
                    top_start: fragment.top_start,
                    bottom_start: fragment.bottom_start,
                });
                if index == request.fragment_index {
                    selected = Some(fragment);
                }
            }
            Ok(geometries)
        },
    )?;
    if !inspection.products_determined {
        return Err(bad("bounded actual digest products are not determined"));
    }
    let fragment =
        selected.ok_or_else(|| bad("selected fragment does not exist in the actual digest"))?;
    let mut reference = fragment.duplex.duplex.molecule(format!(
        "Saved digest fragment {}",
        request.fragment_index + 1
    ))?;
    let geometry = chemistry::DigestFragment {
        duplex: fragment.duplex.duplex.clone(),
        top_start: fragment.top_start,
        bottom_start: fragment.bottom_start,
    };
    let mut findings = vec![];
    verification::project_digest_annotations(source, &geometry, &mut reference, &mut findings)?;
    reference
        .metadata
        .insert(CONTEXT_KEY.into(), context_revision.into());
    reference.metadata.insert(
        catalog_digest::IMPORTED_RECIPE_KEY.into(),
        serde_json::to_string(&request.selection).map_err(|_| Error::TooLarge)?,
    );
    reference.metadata.insert("source:run_digest_history:scope".into(),"Target-free actual cut derivation. Source recipe, chemistry and modifications require native replay; detached metadata is not material authority or independent-target verification.".into());
    if !findings.is_empty() {
        reference.metadata.insert(
            "source:run_digest_history:annotation_findings".into(),
            serde_json::to_string(&findings).map_err(|_| Error::TooLarge)?,
        );
    }
    reference.validate()?;
    fragment.duplex.validate()?;
    Ok(FreshFragment {
        reference,
        tracked: fragment.duplex,
        origin_incomplete: replay.incomplete,
        origin_findings,
        origin_findings_total,
        origin_findings_truncated,
    })
}
pub fn prepare_install(w: &Workspace, request: &Request) -> Result<PreparedInstallation, Error> {
    let context_revision = installation_revision(w, request)?;
    let fingerprint = run_lineage::workspace_fingerprint(w)?;
    let source = run_material::resolve(w, &request.source)?;
    let fresh = fresh(source, request, &context_revision)?;
    let native = restriction_history::prepare_fresh_install(
        w,
        &request.source,
        fresh.reference.clone(),
        context_revision.clone(),
        format!(
            "Target-free saved-run digest fragment {} · originating event {} · catalog-pinned complete-cleavage assumption; no independent fragment target was assessed",
            request.fragment_index + 1,
            request.source.event_id
        ),
    )?;
    let product = EventBinding {
        document_id: native.document_id,
        event_id: native.event_id,
        revision: native.product_revision,
    };
    let record = Record {
        schema_version: 1,
        engine_version: ENGINE_VERSION.into(),
        scientific_engine: verification::ENGINE_VERSION.into(),
        source_workspace_fingerprint: fingerprint,
        context_revision: context_revision.clone(),
        request: request.clone(),
        product: product.clone(),
        product_state_revision: modified_duplex::summary(&fresh.tracked)?.state_revision,
        origin_incomplete: fresh.origin_incomplete,
    };
    record.validate_documentary()?;
    let bytes = encoded(&record, MAX_RECORD_BYTES)?;
    let record_hash = blake3::hash(&bytes).to_hex().to_string();
    let mut workspace = native.workspace;
    if workspace
        .attachments
        .insert(
            record_hash.clone(),
            WorkspaceAttachment {
                kind: WORKSPACE_ATTACHMENT_KIND.into(),
                name: format!("Saved digest fragment {}", request.fragment_index + 1),
                media_type: "application/json".into(),
                bytes,
            },
        )
        .is_some()
    {
        return Err(bad("derivation attachment collision"));
    }
    workspace.validate()?;
    validate_workspace_records(&workspace)?;
    let receipt = Receipt {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        record_hash,
        source: request.source.clone(),
        product,
        origin_incomplete: fresh.origin_incomplete,
        origin_findings: fresh.origin_findings,
        origin_findings_total: fresh.origin_findings_total,
        origin_findings_truncated: fresh.origin_findings_truncated,
        duplex: fresh.tracked.duplex.summary()?,
        modifications: modified_duplex::summary(&fresh.tracked)?,
    };
    Ok(PreparedInstallation { workspace, receipt })
}
impl Receipt {
    /// Bounded documentary binding only; this never reruns the origin or digest.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if self.engine_version != ENGINE_VERSION
            || !hex(expected)
            || self.context_revision != expected
            || !hex(&self.record_hash)
            || self.origin_findings.len() != self.origin_findings_total.min(128)
            || self
                .origin_findings
                .iter()
                .any(|s| s.chars().count() > 2048)
            || (self.origin_findings_total > self.origin_findings.len()
                && !self.origin_findings_truncated)
            || self.origin_findings_total > crate::run_lineage::MAX_RUN_BYTES
            || self.duplex.top_length == 0
            || self.duplex.bottom_length == 0
            || self.duplex.top_length > 1_000_000
            || self.duplex.bottom_length > 1_000_000
            || self.modifications.top_length != self.duplex.top_length
            || self.modifications.bottom_length != self.duplex.bottom_length
            || !hex(&self.modifications.state_revision)
        {
            return Err(bad("invalid installed receipt/context or bounded findings"));
        }
        run_material::event_shape(&self.source)?;
        run_material::event_shape(&self.product)?;
        encoded(self, MAX_RECORD_BYTES)?;
        Ok(())
    }
    pub fn validate_archive(&self, w: &Workspace) -> Result<(), Error> {
        self.validate_binding(&self.context_revision)?;
        w.validate()?;
        validate_workspace_records(w)?;
        let attachment = w
            .attachments
            .get(&self.record_hash)
            .filter(|a| a.kind == WORKSPACE_ATTACHMENT_KIND)
            .ok_or_else(|| bad("installed derivation record unavailable"))?;
        let record = load_record(&attachment.bytes)?;
        if record.context_revision != self.context_revision
            || record.product != self.product
            || record.request.source != self.source
            || record.origin_incomplete != self.origin_incomplete
            || record.product_state_revision != self.modifications.state_revision
        {
            return Err(bad(
                "installed receipt differs from exact derivation record",
            ));
        }
        Ok(())
    }
}
impl Installed {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.receipt.validate_binding(expected)?;
        let w = workspace::load_workspace(&self.archive)?;
        self.receipt.validate_archive(&w)
    }
}
pub fn install(w: &Workspace, request: &Request) -> Result<Installed, Error> {
    let prepared = prepare_install(w, request)?;
    Ok(Installed {
        archive: workspace::save_workspace(&prepared.workspace)?,
        receipt: prepared.receipt,
    })
}
/// Safe structural hook for workspace validators. Never replays science or
/// treats a marker/attachment as an authoritative physical product.
pub fn validate_workspace_records(w: &Workspace) -> Result<(), Error> {
    let mut products = BTreeSet::new();
    let mut count = 0;
    for (hash, a) in w
        .attachments
        .iter()
        .filter(|(_, a)| a.kind == WORKSPACE_ATTACHMENT_KIND)
    {
        count += 1;
        if count > MAX_RECORDS
            || a.media_type != "application/json"
            || blake3::hash(&a.bytes).to_hex().as_str() != hash
        {
            return Err(bad("derivation attachment count/media/hash changed"));
        }
        let record = load_record(&a.bytes)?;
        if !products.insert(record.product.event_id.clone()) {
            return Err(bad("multiple derivations claim the same product event"));
        }
        attest(w, &record.request)?;
        let reference = run_material::resolve(w, &record.product)?;
        if reference.metadata.get(CONTEXT_KEY) != Some(&record.context_revision)
            || !reference.metadata.contains_key("duplex_model")
        {
            return Err(bad(
                "derived reference lacks its restrictive native context",
            ));
        }
        let event = w
            .history
            .iter()
            .find(|e| e.id == record.product.event_id)
            .ok_or_else(|| bad("derived product event missing"))?;
        if event.parents != vec![record.request.source.event_id.clone()] {
            return Err(bad("derived event no longer names its exact source parent"));
        }
    }
    Ok(())
}
fn record_for(w: &Workspace, product: &EventBinding) -> Result<Record, Error> {
    run_material::event_shape(product)?;
    let mut found = None;
    for a in w
        .attachments
        .values()
        .filter(|a| a.kind == WORKSPACE_ATTACHMENT_KIND)
    {
        let r = load_record(&a.bytes)?;
        if &r.product == product {
            if found.is_some() {
                return Err(bad("ambiguous derived product event"));
            }
            found = Some(r);
        }
    }
    found.ok_or_else(|| bad("exact derived product record unavailable"))
}
/// Close the target-free recipe with its exact selected source snapshot. The
/// returned serializable object is still documentary; replay grants authority.
pub fn capsule(w: &Workspace, product: &EventBinding) -> Result<Capsule, Error> {
    let _fresh = reconstruct(w, product)?;
    let record = record_for(w, product)?;
    let source = run_material::resolve(w, &record.request.source)?;
    encoded(
        &(1u32, &record, source),
        MAX_CAPSULE_BYTES.saturating_sub(128),
    )?;
    let source = source.clone();
    let capsule = Capsule {
        schema_version: 1,
        record,
        source,
    };
    capsule.validate_documentary()?;
    Ok(capsule)
}
impl Capsule {
    pub fn validate_documentary(&self) -> Result<(), Error> {
        crate::run_material_capsule::digest_capsule_preflight(self)
    }
}

/// Receiving workspace/native-host check before selecting a serialized capsule.
/// This verifies retained ancestry, not chemical state.
pub fn validate_capsule_for_workspace(
    w: &Workspace,
    capsule: &Capsule,
    product: &EventBinding,
) -> Result<(), Error> {
    if &capsule.record.product != product {
        return Err(bad(
            "digest material requires its exact saved product event",
        ));
    }
    crate::run_material_capsule::validate_digest_capsule_origins(w, capsule)
}
/// Closed recipe replay for a future typed FromDerivedDigest operation. Never
/// trusts a cached duplex, material summary, fragment index chosen by target DNA,
/// or an imported verification result.
pub fn replay(capsule: &Capsule, reference: &Molecule) -> Result<FreshFragment, Error> {
    capsule.validate_documentary()?;
    crate::model::bounds::model_preflight(&[reference])?;
    replay_fields(&capsule.record, &capsule.source, reference)
}
fn replay_fields(
    record: &Record,
    source: &Molecule,
    reference: &Molecule,
) -> Result<FreshFragment, Error> {
    if record.scientific_engine != verification::ENGINE_VERSION {
        return Err(bad(
            "historical scientific engine remains readable evidence but cannot grant current digest material",
        ));
    }
    let fresh = fresh(source, &record.request, &record.context_revision)?;
    if modified_duplex::summary(fresh.tracked())?.state_revision != record.product_state_revision
        || fresh.origin_incomplete != record.origin_incomplete
        || reference.revision() != record.product.revision
        || reference != fresh.reference()
    {
        return Err(bad(
            "derived product/current origin assessment differs from fresh actual cleavage",
        ));
    }
    Ok(fresh)
}
pub fn reconstruct(w: &Workspace, product: &EventBinding) -> Result<FreshFragment, Error> {
    w.validate()?;
    validate_workspace_records(w)?;
    let record = record_for(w, product)?;
    require_current_origins(w, &record.request)?;
    let source = run_material::resolve(w, &record.request.source)?;
    let reference = run_material::resolve(w, product)?;
    replay_fields(&record, source, reference)
}
