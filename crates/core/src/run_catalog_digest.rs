//! Catalog-pinned virtual cleavage of freshly reconstructed saved material.
//! Cleavage remains the existing declared complete-cleavage model, never a
//! methylation permission assessment or an experimental reaction claim.
use crate::{
    Error, Molecule, assembly_recipe,
    catalog_digest::{self, CatalogRestrictionInspection, CatalogSelection, catalog},
    chemistry,
    modified_duplex::{self, ModificationSummary},
    restriction_inspection::{self, InspectionOptions},
    run_material::Capsule,
    verification::{self, Status},
};
use serde::{Deserialize, Serialize};
pub const ENGINE_VERSION: &str = "run-catalog-digest-1";
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Inspection{pub engine_version:String,pub context_revision:String,pub origin_revision:String,pub origin_status:Status,pub source_duplex:chemistry::DuplexSummary,pub source_modifications:ModificationSummary,pub fragment_modifications:Vec<ModificationSummary>,pub inspection:CatalogRestrictionInspection,pub findings:Vec<String>}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("saved material catalog digest: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn context(origin: &str, batch: &str) -> Result<String, Error> {
    assembly_recipe::hash(&(ENGINE_VERSION, verification::ENGINE_VERSION, origin, batch))
}
pub fn source_preflight(
    c: &catalog::Catalog,
    s: &CatalogSelection,
    role: &str,
    m: &Molecule,
    origin: &Capsule,
    o: &InspectionOptions,
) -> Result<(), Error> {
    origin.validate_documentary()?;
    catalog_digest::source_preflight(c, s, role, m, o)?;
    if m.revision() != origin.selection.event.revision {
        return Err(bad("saved source differs from the exact originating event"));
    }
    Ok(())
}
pub fn inspection_revision(
    c: &catalog::Catalog,
    s: &CatalogSelection,
    role: &str,
    m: &Molecule,
    origin: &Capsule,
    o: &InspectionOptions,
) -> Result<String, Error> {
    source_preflight(c, s, role, m, origin, o)?;
    context(
        &assembly_recipe::hash(origin)?,
        &catalog_digest::batch_revision(c, s, role, m, o)?,
    )
}
pub fn inspect(
    c: &catalog::Catalog,
    s: &CatalogSelection,
    role: &str,
    m: &Molecule,
    origin: &Capsule,
    o: &InspectionOptions,
) -> Result<Inspection, Error> {
    let context_revision = inspection_revision(c, s, role, m, origin, o)?;
    let fresh = crate::run_material_capsule::replay(origin, m)?;
    let origin_report = fresh.evidence.decode_report()?;
    let source_modifications = modified_duplex::summary(&fresh.tracked)?;
    let source_duplex = fresh.tracked.duplex.summary()?;
    let enzymes = s.resolved_enzymes(c)?;
    let mut fragment_modifications = vec![];
    let inner = restriction_inspection::inspect_with_digest(role, m, &enzymes, o, || {
        let fragments = modified_duplex::digest(&fresh.tracked, &enzymes)?;
        let mut result = Vec::with_capacity(fragments.len());
        for fragment in fragments {
            fragment_modifications.push(modified_duplex::summary(&fragment.duplex)?);
            result.push(chemistry::DigestFragment {
                duplex: fragment.duplex.duplex,
                top_start: fragment.top_start,
                bottom_start: fragment.bottom_start,
            });
        }
        Ok(result)
    })?;
    if !inner.products_determined {
        fragment_modifications.clear();
    }
    let inspection = CatalogRestrictionInspection {
        engine_version: catalog_digest::ENGINE_VERSION.into(),
        selection: s.clone(),
        batch_revision: catalog_digest::batch_revision(c, s, role, m, o)?,
        inspection: inner,
    };
    let mut findings=vec!["Recognition and products use the actual freshly replayed saved strand material. Complete cleavage is an explicit virtual-digest assumption, not measured enzyme activity or methylation permission. Positional modifications, generation declarations and terminal groups are retained; imported summaries never establish material authority. Existing nick propagation remains unsupported and withholds products while recognition sites remain inspectable.".into()];
    if fresh.incomplete {
        findings.push("The originating run remains incomplete. Its retained assumptions and unknown/unapplied conditions are not resolved by displaying a selected digest fragment.".into());
    }
    if source_modifications.known_modified > 0 || source_modifications.unknown_eligible > 0 {
        findings.push("Inherited modified or chemically unknown eligible bases are present. This catalog complete-cleavage profile does not assess methylation susceptibility; product geometry is conditional on the declared cleavage assumption. Use the existing positional checked-digest profile where applicable before asserting cleavage permission.".into());
    }
    findings.extend(
        origin_report
            .findings
            .into_iter()
            .map(|s| format!("Origin finding: {s}")),
    );
    findings.extend(
        origin_report
            .assumptions
            .into_iter()
            .map(|s| format!("Origin assumption: {s}")),
    );
    let report = Inspection {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        origin_revision: assembly_recipe::hash(origin)?,
        origin_status: origin_report.status,
        source_duplex,
        source_modifications,
        fragment_modifications,
        inspection,
        findings,
    };
    report.validate_binding(&report.context_revision)?;
    Ok(report)
}
fn summary_shape(s: &ModificationSummary) -> Result<(), Error> {
    let total = s
        .top_length
        .checked_add(s.bottom_length)
        .ok_or(Error::TooLarge)?;
    let calls = s
        .known_modified
        .checked_add(s.known_unmodified)
        .and_then(|n| n.checked_add(s.unknown_eligible))
        .ok_or(Error::TooLarge)?;
    if s.engine_version != modified_duplex::ENGINE_VERSION
        || !hex(&s.state_revision)
        || s.top_length > crate::model::MAX_SEQUENCE_LEN
        || s.bottom_length > crate::model::MAX_SEQUENCE_LEN
        || calls > total
        || s.generation.len() > 1000
        || s.generation.iter().any(|g| !hex(&g.source_revision))
    {
        return Err(bad("invalid bounded modification summary"));
    }
    Ok(())
}
impl Inspection {
    /// Readonly documentary consistency. This does not replay material or assess cleavage.
    pub fn validate(&self) -> Result<(), Error> {
        self.inspection.validate()?;
        summary_shape(&self.source_modifications)?;
        crate::digest_cut_selection::source_geometry_preflight(
            &self.inspection.inspection,
            &self.source_duplex,
        )?;
        if self.engine_version != ENGINE_VERSION
            || !hex(&self.origin_revision)
            || self.origin_status == Status::Invalid
            || self.context_revision
                != context(&self.origin_revision, &self.inspection.batch_revision)?
            || self.source_modifications.top_length != self.inspection.inspection.source_length
            || self.source_modifications.bottom_length != self.source_duplex.bottom_length
            || self.fragment_modifications.len() != self.inspection.inspection.fragments.len()
            || self.findings.len() > 20_000
            || self.findings.iter().any(|s| s.len() > 131_072)
        {
            return Err(bad(
                "inspection context, source or fragment correspondence changed",
            ));
        }
        for (summary, fragment) in self
            .fragment_modifications
            .iter()
            .zip(&self.inspection.inspection.fragments)
        {
            summary_shape(summary)?;
            if summary.top_length != fragment.duplex.top_length
                || summary.bottom_length != fragment.duplex.bottom_length
                || summary.generation != self.source_modifications.generation
            {
                return Err(bad(
                    "fragment modification/source geometry correspondence changed",
                ));
            }
        }
        assembly_recipe::hash(self)?;
        Ok(())
    }
    pub fn selected_fragment(
        &self,
        choice: &crate::digest_cut_selection::Choice,
    ) -> Result<&crate::restriction_inspection::FragmentSummary, Error> {
        self.validate()?;
        crate::digest_cut_selection::selected_fragment_with_source(
            &self.inspection,
            choice,
            &self.source_duplex,
        )
    }
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.validate()?;
        if self.context_revision != expected {
            return Err(bad("stale origin/catalog/source/options context"));
        }
        Ok(())
    }
    pub fn validate_bindings_for(
        &self,
        c: &catalog::Catalog,
        s: &CatalogSelection,
        role: &str,
        m: &Molecule,
        origin: &Capsule,
        o: &InspectionOptions,
    ) -> Result<(), Error> {
        self.validate_binding(&inspection_revision(c, s, role, m, origin, o)?)?;
        self.inspection.validate_bindings_for(c, s, role, m, o)?;
        if self.origin_revision != assembly_recipe::hash(origin)? {
            return Err(bad("saved origin declaration changed"));
        }
        Ok(())
    }
}
