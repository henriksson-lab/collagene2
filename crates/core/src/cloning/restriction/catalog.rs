//! Catalog-pinned recognition/assumed-cleavage workflows using the existing scientific backend.
//! Identity pins and supplier observations are documentary producer data, never reaction validation.
use crate::cloning::restriction::inspection::InspectionOptions;
use crate::{
    cloning::restriction::{
        history::{self as restriction_history, InstalledDigest, PreparedInstallation},
        inspection::{
            self as restriction_inspection, DigestProduct, EnzymeSpec, RestrictionInspection,
        },
    },
    documents::workspace::{self, Workspace},
    material::lineage::{self as run_lineage, EventBinding as SourceBinding},
    sequence::analysis::Enzyme,
};
pub use collagene_enzyme_catalog as catalog;
use collagene_model::{Error, Molecule};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeSet, io::Write};
pub const ENGINE_VERSION: &str = "catalog-digest-1";
pub const SCHEMA_VERSION: u16 = 1;
pub const MAX_SELECTION_BYTES: usize = 32 * 1024;
pub const RECIPE_KEY: &str = "restriction_catalog:selection";
pub const IMPORTED_RECIPE_KEY: &str = "source:restriction_catalog:selection";
pub const ENGINE_KEY: &str = "restriction_catalog:engine";
pub const PRODUCT_REVISION_KEY: &str = "restriction_catalog:batch_revision";
macro_rules! dto{($($x:tt)*)=>{#[derive(Debug,Clone,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($x)*};}
dto! {#[derive(PartialEq,Eq)] pub struct SelectedEnzyme{pub id:String,pub enzyme:EnzymeSpec}}
dto! {#[derive(PartialEq,Eq)] pub struct CatalogSelection{pub schema_version:u16,pub catalog_edition:String,pub catalog_content_hash:String,#[serde(deserialize_with="collagene_model::decode::bounded_vec::<_, SelectedEnzyme, 64>")]pub enzymes:Vec<SelectedEnzyme>}}
dto! {pub struct CatalogRestrictionInspection{pub engine_version:String,pub selection:CatalogSelection,pub batch_revision:String,pub inspection:RestrictionInspection}}
dto! {pub struct CatalogDigestProduct{pub engine_version:String,pub selection:CatalogSelection,pub batch_revision:String,pub product:DigestProduct}}
// batch_revision here is the installation request revision, including the original
// workspace fingerprint. The molecule's PRODUCT_REVISION_KEY is the separate product revision.
dto! {pub struct CatalogDigestInstalled{pub engine_version:String,pub selection:CatalogSelection,pub batch_revision:String,pub installed:InstalledDigest}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("catalog digest: {s}"))
}
fn hash<T: Serialize>(v: &T) -> Result<String, Error> {
    Ok(
        blake3::hash(&serde_json::to_vec(v).map_err(|_| Error::TooLarge)?)
            .to_hex()
            .to_string(),
    )
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || matches!(b, b'a'..=b'f'))
}
fn text(s: &str, n: usize) -> bool {
    s.len() <= n && !s.trim().is_empty() && !s.chars().any(char::is_control)
}
fn engine(s: &str) -> bool {
    text(s, 128) && s.is_ascii()
}
struct Bound(usize);
impl Write for Bound {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.0 = self
            .0
            .checked_sub(b.len())
            .ok_or_else(|| std::io::Error::other("bounded catalog selection"))?;
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn selection_size(s: &CatalogSelection) -> Result<(), Error> {
    serde_json::to_writer(Bound(MAX_SELECTION_BYTES), s).map_err(|_| Error::TooLarge)
}
fn spec(e: &catalog::Enzyme) -> EnzymeSpec {
    EnzymeSpec {
        name: e.name.clone(),
        recognition: e.geometry.recognition.clone(),
        top_cut: e.geometry.top_cut,
        bottom_cut: e.geometry.bottom_cut,
    }
}
impl CatalogSelection {
    /// Historical shape only. Membership and supplier authenticity are not established here.
    pub fn validate(&self) -> Result<(), Error> {
        if self.schema_version != SCHEMA_VERSION
            || !text(&self.catalog_edition, 128)
            || !hex(&self.catalog_content_hash)
            || self.enzymes.is_empty()
            || self.enzymes.len() > 64
        {
            return Err(bad("invalid catalog pin or selection capacity"));
        }
        let mut ids = BTreeSet::new();
        let mut names = BTreeSet::new();
        for selected in &self.enzymes {
            let e = &selected.enzyme;
            if !text(&selected.id, 96)
                || !ids.insert(selected.id.as_str())
                || !text(&e.name, 64)
                || !names.insert(e.name.to_ascii_lowercase())
                || e.recognition.is_empty()
                || e.recognition.len() > 64
                || !e
                    .recognition
                    .bytes()
                    .all(|b| b"ACGTRYSWKMBDHVN".contains(&b))
                || !(-1024..=1024).contains(&e.top_cut)
                || !(-1024..=1024).contains(&e.bottom_cut)
            {
                return Err(bad(
                    "invalid, duplicate or unsupported selected identity/geometry",
                ));
            }
        }
        selection_size(self)
    }
    /// Exact full catalog pin and ordered ID/geometry/name membership; never alias expansion.
    pub fn validate_for(&self, c: &catalog::Catalog) -> Result<(), Error> {
        self.validate()?;
        if self.catalog_edition != c.edition().edition
            || self.catalog_content_hash != c.content_hash()
        {
            return Err(Error::Format(
                "catalog digest: catalog edition/content changed; select again explicitly".into(),
            ));
        }
        for selected in &self.enzymes {
            let record = c
                .edition()
                .enzymes
                .iter()
                .find(|e| e.id == selected.id)
                .ok_or_else(|| {
                    Error::Format("catalog digest: selected exact enzyme ID is absent".into())
                })?;
            if selected.enzyme != spec(record) {
                return Err(Error::Format(
                    "catalog digest: selected identity/name/geometry differs from the pinned catalog".into(),
                ));
            }
        }
        Ok(())
    }
    pub fn resolved_enzymes(&self, c: &catalog::Catalog) -> Result<Vec<Enzyme>, Error> {
        self.validate_for(c)?;
        Ok(self
            .enzymes
            .iter()
            .map(|s| Enzyme {
                name: s.enzyme.name.clone(),
                recognition: s.enzyme.recognition.clone(),
                top_cut: s.enzyme.top_cut,
                bottom_cut: s.enzyme.bottom_cut,
            })
            .collect())
    }
}
pub fn builtin_catalog() -> Result<catalog::Catalog, Error> {
    catalog::Catalog::builtin().map_err(|e| bad(&e.to_string()))
}
/// Exact canonical enzyme IDs only. Names, aliases, vendor codes and formulations
/// must be resolved by an explicit user choice; no formulation conditions imply cleavage.
pub fn select(c: &catalog::Catalog, ids: &[String]) -> Result<CatalogSelection, Error> {
    if ids.is_empty() || ids.len() > 64 || ids.iter().any(|id| !text(id, 96)) {
        return Err(bad("select 1–64 exact enzyme IDs"));
    }
    let mut enzymes = Vec::with_capacity(ids.len());
    for id in ids {
        let record = c
            .edition()
            .enzymes
            .iter()
            .find(|e| &e.id == id)
            .ok_or_else(|| {
                bad("exact enzyme ID required; aliases/formulations are not implicit selections")
            })?;
        enzymes.push(SelectedEnzyme {
            id: id.clone(),
            enzyme: spec(record),
        });
    }
    let selected = CatalogSelection {
        schema_version: SCHEMA_VERSION,
        catalog_edition: c.edition().edition.clone(),
        catalog_content_hash: c.content_hash().into(),
        enzymes,
    };
    selected.validate_for(c)?;
    Ok(selected)
}
pub fn source_preflight(
    c: &catalog::Catalog,
    s: &CatalogSelection,
    role: &str,
    m: &Molecule,
    o: &InspectionOptions,
) -> Result<(), Error> {
    s.validate_for(c)?;
    restriction_inspection::source_preflight(role, m, &s.resolved_enzymes(c)?, o)
}
fn inspection_hash(engine: &str, s: &CatalogSelection, base: &str) -> Result<String, Error> {
    hash(&("catalog-restriction-inspection-1", engine, s, base))
}
fn product_hash(engine: &str, s: &CatalogSelection, base: &str) -> Result<String, Error> {
    hash(&("catalog-restriction-product-1", engine, s, base))
}
pub fn batch_revision(
    c: &catalog::Catalog,
    s: &CatalogSelection,
    role: &str,
    m: &Molecule,
    o: &InspectionOptions,
) -> Result<String, Error> {
    source_preflight(c, s, role, m, o)?;
    inspection_hash(
        ENGINE_VERSION,
        s,
        &restriction_inspection::batch_revision(role, m, &s.resolved_enzymes(c)?, o)?,
    )
}
pub fn product_revision(
    c: &catalog::Catalog,
    s: &CatalogSelection,
    role: &str,
    m: &Molecule,
    o: &InspectionOptions,
    index: usize,
) -> Result<String, Error> {
    source_preflight(c, s, role, m, o)?;
    product_hash(
        ENGINE_VERSION,
        s,
        &restriction_inspection::product_revision(role, m, &s.resolved_enzymes(c)?, o, index)?,
    )
}
pub fn inspect(
    c: &catalog::Catalog,
    s: &CatalogSelection,
    role: &str,
    m: &Molecule,
    o: &InspectionOptions,
) -> Result<CatalogRestrictionInspection, Error> {
    source_preflight(c, s, role, m, o)?;
    let inspection = restriction_inspection::inspect(role, m, &s.resolved_enzymes(c)?, o)?;
    let r = CatalogRestrictionInspection {
        engine_version: ENGINE_VERSION.into(),
        selection: s.clone(),
        batch_revision: inspection_hash(ENGINE_VERSION, s, &inspection.batch_revision)?,
        inspection,
    };
    r.validate()?;
    Ok(r)
}
pub fn product(
    c: &catalog::Catalog,
    s: &CatalogSelection,
    role: &str,
    m: &Molecule,
    o: &InspectionOptions,
    index: usize,
) -> Result<CatalogDigestProduct, Error> {
    source_preflight(c, s, role, m, o)?;
    let mut product = restriction_inspection::product(role, m, &s.resolved_enzymes(c)?, o, index)?;
    let revision = product_hash(ENGINE_VERSION, s, &product.batch_revision)?;
    product.molecule.metadata.insert(
        RECIPE_KEY.into(),
        serde_json::to_string(s).map_err(|_| Error::TooLarge)?,
    );
    product
        .molecule
        .metadata
        .insert(ENGINE_KEY.into(), ENGINE_VERSION.into());
    product
        .molecule
        .metadata
        .insert(PRODUCT_REVISION_KEY.into(), revision.clone());
    let p = CatalogDigestProduct {
        engine_version: ENGINE_VERSION.into(),
        selection: s.clone(),
        batch_revision: revision,
        product,
    };
    p.validate()?;
    crate::protocol::save_document(&p.product.molecule)?;
    Ok(p)
}
/// Read bounded unverified producer data. This cannot establish installed catalog
/// membership, source chemistry, cleavage, or current verification.
pub fn unverified_selection(m: &Molecule) -> Result<Option<CatalogSelection>, Error> {
    let Some(raw) = m
        .metadata
        .get(RECIPE_KEY)
        .or_else(|| m.metadata.get(IMPORTED_RECIPE_KEY))
    else {
        return Ok(None);
    };
    if raw.len() > MAX_SELECTION_BYTES {
        return Err(Error::TooLarge);
    }
    let s: CatalogSelection =
        serde_json::from_str(raw).map_err(|_| bad("invalid catalog selection JSON"))?;
    s.validate()?;
    Ok(Some(s))
}
fn metadata(m: &Molecule, s: &CatalogSelection, engine: &str, revision: &str) -> Result<(), Error> {
    if !m.metadata.contains_key(RECIPE_KEY)
        || unverified_selection(m)?.as_ref() != Some(s)
        || m.metadata.get(ENGINE_KEY).map(String::as_str) != Some(engine)
        || m.metadata.get(PRODUCT_REVISION_KEY).map(String::as_str) != Some(revision)
    {
        return Err(bad("producer catalog metadata differs from request/recipe"));
    }
    Ok(())
}
impl CatalogRestrictionInspection {
    pub fn validate(&self) -> Result<(), Error> {
        self.selection.validate()?;
        self.inspection.validate()?;
        if !engine(&self.engine_version)
            || self.inspection.enzymes
                != self
                    .selection
                    .enzymes
                    .iter()
                    .map(|e| e.enzyme.clone())
                    .collect::<Vec<_>>()
            || self.batch_revision
                != inspection_hash(
                    &self.engine_version,
                    &self.selection,
                    &self.inspection.batch_revision,
                )?
        {
            return Err(bad("inspection catalog binding mismatch"));
        }
        Ok(())
    }
    pub fn validate_bindings_for(
        &self,
        c: &catalog::Catalog,
        s: &CatalogSelection,
        role: &str,
        m: &Molecule,
        o: &InspectionOptions,
    ) -> Result<(), Error> {
        source_preflight(c, s, role, m, o)?;
        self.validate()?;
        if self.engine_version != ENGINE_VERSION
            || &self.selection != s
            || self.batch_revision != batch_revision(c, s, role, m, o)?
        {
            return Err(bad("stale inspection catalog/source/options"));
        }
        self.inspection
            .validate_bindings_for(role, m, &s.resolved_enzymes(c)?, o)
    }
}
impl CatalogDigestProduct {
    pub fn validate(&self) -> Result<(), Error> {
        self.selection.validate()?;
        self.product.validate()?;
        if !engine(&self.engine_version)
            || self.product.enzymes
                != self
                    .selection
                    .enzymes
                    .iter()
                    .map(|e| e.enzyme.clone())
                    .collect::<Vec<_>>()
            || self.batch_revision
                != product_hash(
                    &self.engine_version,
                    &self.selection,
                    &self.product.batch_revision,
                )?
        {
            return Err(bad("product catalog binding mismatch"));
        }
        metadata(
            &self.product.molecule,
            &self.selection,
            &self.engine_version,
            &self.batch_revision,
        )
    }
    pub fn validate_bindings_for(
        &self,
        c: &catalog::Catalog,
        s: &CatalogSelection,
        role: &str,
        m: &Molecule,
        o: &InspectionOptions,
        index: usize,
    ) -> Result<(), Error> {
        source_preflight(c, s, role, m, o)?;
        self.validate()?;
        if self.engine_version != ENGINE_VERSION
            || &self.selection != s
            || self.batch_revision != product_revision(c, s, role, m, o, index)?
        {
            return Err(bad("stale product catalog/source/options/index"));
        }
        self.product
            .validate_bindings_for(role, m, &s.resolved_enzymes(c)?, o, index)
    }
}
pub fn installation_revision(
    c: &catalog::Catalog,
    s: &CatalogSelection,
    w: &Workspace,
    b: &SourceBinding,
    o: &InspectionOptions,
    index: usize,
) -> Result<String, Error> {
    s.validate_for(c)?;
    let m = restriction_history::current_source(w, b)?;
    let revision = product_revision(c, s, &restriction_history::source_role(b)?, m, o, index)?;
    hash(&(
        "catalog-restriction-installation-1",
        ENGINE_VERSION,
        run_lineage::workspace_fingerprint(w)?,
        b,
        revision,
    ))
}
pub fn prepare_install(
    c: &catalog::Catalog,
    s: &CatalogSelection,
    w: &Workspace,
    b: &SourceBinding,
    o: &InspectionOptions,
    index: usize,
) -> Result<PreparedInstallation, Error> {
    let revision = installation_revision(c, s, w, b, o, index)?;
    let m = restriction_history::current_source(w, b)?;
    let p = product(c, s, &restriction_history::source_role(b)?, m, o, index)?;
    p.validate_bindings_for(c, s, &restriction_history::source_role(b)?, m, o, index)?;
    restriction_history::prepare_fresh_install(
        w,
        b,
        p.product.molecule,
        revision,
        format!(
            "Catalog virtual digest fragment {} · catalog {} / {} · source event {}; complete-cleavage assumptions only; catalog evidence does not validate reaction conditions",
            index + 1,
            s.catalog_edition,
            s.catalog_content_hash,
            b.event_id
        ),
    )
}
pub fn install(
    c: &catalog::Catalog,
    s: &CatalogSelection,
    w: &Workspace,
    b: &SourceBinding,
    o: &InspectionOptions,
    index: usize,
) -> Result<CatalogDigestInstalled, Error> {
    let p = prepare_install(c, s, w, b, o, index)?;
    let installed = InstalledDigest {
        archive: workspace::save_workspace(&p.workspace)?,
        source: p.source,
        document_id: p.document_id,
        event_id: p.event_id,
        product_revision: p.product_revision,
        batch_revision: p.batch_revision.clone(),
    };
    let result = CatalogDigestInstalled {
        engine_version: ENGINE_VERSION.into(),
        selection: s.clone(),
        batch_revision: p.batch_revision,
        installed,
    };
    result.validate_archive(&p.workspace)?;
    Ok(result)
}
impl CatalogDigestInstalled {
    /// Envelope validation only; catalog imports and producer metadata confer no scientific authority.
    pub fn validate(&self) -> Result<(), Error> {
        self.selection.validate()?;
        self.installed.validate()?;
        if !engine(&self.engine_version)
            || !hex(&self.batch_revision)
            || self.batch_revision != self.installed.batch_revision
        {
            return Err(bad("invalid installation catalog envelope"));
        }
        Ok(())
    }
    /// Cheap binding check for a decoded archive returned by a fresh private worker.
    /// This does not establish science for externally supplied archive/recipe claims.
    pub fn validate_archive(&self, w: &Workspace) -> Result<(), Error> {
        self.validate()?;
        self.installed.validate_archive(w)?;
        let m = w
            .revisions
            .get(&self.installed.product_revision)
            .ok_or_else(|| bad("catalog product revision absent"))?;
        let recipe = restriction_inspection::unverified_recipe(m)?
            .ok_or_else(|| bad("catalog product recipe absent"))?;
        let selected: Vec<_> = self
            .selection
            .enzymes
            .iter()
            .map(|e| e.enzyme.clone())
            .collect();
        if recipe.source_role != restriction_history::source_role(&self.installed.source)?
            || recipe.source_revision != self.installed.source.revision
            || recipe.enzymes != selected
        {
            return Err(bad(
                "catalog producer recipe has another source or selected geometry",
            ));
        }
        metadata(
            m,
            &self.selection,
            &self.engine_version,
            &product_hash(
                &self.engine_version,
                &self.selection,
                &recipe.batch_revision,
            )?,
        )
    }
    pub fn validate_bindings_for(
        &self,
        c: &catalog::Catalog,
        s: &CatalogSelection,
        w: &Workspace,
        b: &SourceBinding,
        o: &InspectionOptions,
        index: usize,
    ) -> Result<(), Error> {
        s.validate_for(c)?;
        self.validate()?;
        if self.engine_version != ENGINE_VERSION
            || &self.selection != s
            || &self.installed.source != b
            || self.batch_revision != installation_revision(c, s, w, b, o, index)?
        {
            return Err(bad("stale catalog installation request binding"));
        }
        Ok(())
    }
    /// Cheap current-request versus decoded-recipe binding, without site matching
    /// or chemistry. Use only for a fresh private worker response plus the host's
    /// current workspace/epoch guard; arbitrary imported claims remain unverified.
    pub fn validate_archive_bindings_for(
        &self,
        c: &catalog::Catalog,
        s: &CatalogSelection,
        original: &Workspace,
        b: &SourceBinding,
        o: &InspectionOptions,
        index: usize,
        decoded: &Workspace,
    ) -> Result<(), Error> {
        self.validate_bindings_for(c, s, original, b, o, index)?;
        self.validate_archive(decoded)?;
        let m = decoded
            .revisions
            .get(&self.installed.product_revision)
            .ok_or_else(|| bad("decoded product absent"))?;
        let recipe = restriction_inspection::unverified_recipe(m)?
            .ok_or_else(|| bad("decoded recipe absent"))?;
        let source = restriction_history::current_source(original, b)?;
        let role = restriction_history::source_role(b)?;
        let es = s.resolved_enzymes(c)?;
        if recipe.engine_version != restriction_inspection::ENGINE_VERSION
            || recipe.options != *o
            || recipe.fragment_index != index
            || recipe.inspection_revision
                != restriction_inspection::batch_revision(&role, source, &es, o)?
            || recipe.batch_revision
                != restriction_inspection::product_revision(&role, source, &es, o, index)?
        {
            return Err(bad(
                "decoded recipe differs from exact current options/fragment request",
            ));
        }
        metadata(
            m,
            s,
            ENGINE_VERSION,
            &product_revision(c, s, &role, source, o, index)?,
        )
    }
    /// FULL SCIENTIFIC REPLAY for native audit. Do not run on the browser main
    /// thread: reconstruct the pinned product/metadata/event before comparing the
    /// complete resulting workspace. Ordinary UI commits require a fresh worker.
    pub fn validate_archive_for(
        &self,
        c: &catalog::Catalog,
        s: &CatalogSelection,
        original: &Workspace,
        b: &SourceBinding,
        o: &InspectionOptions,
        index: usize,
        decoded: &Workspace,
    ) -> Result<(), Error> {
        self.validate_archive_bindings_for(c, s, original, b, o, index, decoded)?;
        let expected = prepare_install(c, s, original, b, o, index)?;
        if expected.document_id != self.installed.document_id
            || expected.event_id != self.installed.event_id
            || expected.product_revision != self.installed.product_revision
            || run_lineage::workspace_fingerprint(&expected.workspace)?
                != run_lineage::workspace_fingerprint(decoded)?
        {
            return Err(bad(
                "decoded catalog product/workspace differs from fresh scientific replay",
            ));
        }
        Ok(())
    }
}
