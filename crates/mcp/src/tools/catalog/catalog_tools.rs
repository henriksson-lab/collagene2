//! Documentary catalog lookup. Resolution preserves uncertainty and assay context.
use collagene_core::cloning::restriction::catalog as catalog_digest;
use collagene_enzyme_catalog::{self as catalog, Catalog, Lookup, Resolution};
use rmcp::ErrorData;
use schemars::JsonSchema;
use serde::Deserialize;
use serde_json::{Value, json};

#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct CatalogIdentityInput {
    pub catalog_edition: String,
    pub catalog_content_hash: String,
}
#[derive(Debug, Clone, Copy, Default, Deserialize, JsonSchema)]
pub enum LookupMode {
    #[default]
    Search,
    Enzyme,
    Formulation,
    Host,
}
fn default_page() -> usize {
    20
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct LookupQuery {
    /// Search names/aliases or look up an exact enzyme/formulation/host identity. At most 128 bytes.
    pub query: String,
    #[serde(default)]
    pub mode: LookupMode,
    #[serde(default)]
    pub offset: usize,
    #[serde(default = "default_page")]
    pub limit: usize,
    pub catalog_pin: Option<CatalogIdentityInput>,
}
#[derive(Debug, Clone, Copy, Deserialize, JsonSchema)]
#[serde(rename_all = "snake_case")]
pub enum MethylationSystem {
    Dam,
    Dcm,
    CpG,
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum EvidenceCondition {
    IncubationTemperature,
    BufferActivity {
        buffer: String,
        temperature_milli_c: i32,
    },
    Methylation {
        system: MethylationSystem,
    },
}
#[derive(Debug, Clone, Deserialize, JsonSchema)]
#[serde(deny_unknown_fields)]
pub struct EvidenceQuery {
    /// Exact formulation ID/name/alias; enzyme identity alone never supplies formulation facts.
    pub formulation: String,
    pub condition: EvidenceCondition,
    /// Exact assay context, if selected; omitted context can return ContextRequired.
    pub context: Option<String>,
    pub catalog_pin: Option<CatalogIdentityInput>,
}
fn bad(s: impl ToString) -> ErrorData {
    ErrorData::invalid_params(s.to_string(), None)
}
fn bounded(s: &str, n: usize) -> bool {
    s.len() <= n && !s.chars().any(char::is_control)
}
pub fn builtin() -> Result<Catalog, ErrorData> {
    catalog_digest::builtin_catalog().map_err(bad)
}
pub fn resource_uri(c: &Catalog) -> String {
    format!(
        "collagene://enzyme-catalog/v{}/{}/{}",
        c.edition().schema_version,
        c.edition().edition,
        c.content_hash()
    )
}
pub fn identity(c: &Catalog) -> Value {
    json!({"catalog_schema_version":c.edition().schema_version,"catalog_edition":c.edition().edition,"catalog_content_hash":c.content_hash(),"catalog_resource_uri":resource_uri(c)})
}
pub fn resource(c: &Catalog) -> Value {
    json!({"identity":identity(c),"catalog":c.edition(),"authority":"Versioned documentary geometry and formulation/assay evidence, with unresolved facts and citations preserved; supplier references are not measured activity or reaction verification."})
}
fn check_pin(c: &Catalog, pin: Option<&CatalogIdentityInput>) -> Result<(), ErrorData> {
    if pin.is_some_and(|p| {
        p.catalog_edition != c.edition().edition || p.catalog_content_hash != c.content_hash()
    }) {
        return Err(bad(
            "catalog edition/content hash changed; inspect the exact current versioned resource",
        ));
    }
    Ok(())
}
pub fn lookup(args: LookupQuery) -> Result<Value, ErrorData> {
    let c = builtin()?;
    lookup_for(&c, args)
}
fn lookup_for(c: &Catalog, args: LookupQuery) -> Result<Value, ErrorData> {
    check_pin(c, args.catalog_pin.as_ref())?;
    if !bounded(&args.query, 128) || !(1..=100).contains(&args.limit) {
        return Err(bad(
            "catalog query is at most 128 bytes and page limit is 1..100",
        ));
    }
    let mut out = identity(c);
    let body = match args.mode {
        LookupMode::Search => {
            let records = c.search(&args.query, 1000).map_err(bad)?;
            if args.offset > records.len() {
                return Err(bad("catalog search offset exceeds count"));
            }
            let end = args.offset.saturating_add(args.limit).min(records.len());
            json!({"mode":"Search","query":args.query,"count":records.len(),"offset":args.offset,"next_offset":end,"done":end==records.len(),"enzymes":&records[args.offset..end]})
        }
        LookupMode::Enzyme => match c.lookup_enzyme(&args.query) {
            Lookup::Missing => {
                json!({"mode":"Enzyme","query":args.query,"status":"Missing","count":0,"enzymes":[]})
            }
            Lookup::Unique(e) => {
                if args.offset != 0 {
                    return Err(bad("exact unique lookup requires offset 0"));
                }
                json!({"mode":"Enzyme","query":args.query,"status":"Unique","count":1,"enzymes":[e]})
            }
            Lookup::Ambiguous(es) => {
                if args.offset > es.len() {
                    return Err(bad("catalog lookup offset exceeds count"));
                }
                let end = args.offset.saturating_add(args.limit).min(es.len());
                json!({"mode":"Enzyme","query":args.query,"status":"Ambiguous","count":es.len(),"offset":args.offset,"next_offset":end,"done":end==es.len(),"enzymes":&es[args.offset..end]})
            }
        },
        LookupMode::Host => {
            let records = c
                .hosts()
                .iter()
                .filter(|host| args.query.is_empty() || host.id == args.query)
                .collect::<Vec<_>>();
            if args.offset > records.len() {
                return Err(bad("host lookup offset exceeds count"));
            }
            let end = args.offset.saturating_add(args.limit).min(records.len());
            json!({"mode":"Host","query":args.query,"count":records.len(),
                "offset":args.offset,"next_offset":end,"done":end==records.len(),
                "hosts":&records[args.offset..end],
                "authority":"Pinned strain/product-specific expectations only. They do not set Dam/Dcm/CpG chemical calls or establish measured modifications."})
        }
        LookupMode::Formulation => match c.lookup_formulation(&args.query) {
            Lookup::Missing => {
                json!({"mode":"Formulation","query":args.query,"status":"Missing","count":0,"formulations":[]})
            }
            Lookup::Unique(f) => {
                if args.offset != 0 {
                    return Err(bad("exact unique lookup requires offset 0"));
                }
                json!({"mode":"Formulation","query":args.query,"status":"Unique","count":1,"formulations":[f]})
            }
            Lookup::Ambiguous(fs) => {
                if args.offset > fs.len() {
                    return Err(bad("catalog lookup offset exceeds count"));
                }
                let end = args.offset.saturating_add(args.limit).min(fs.len());
                json!({"mode":"Formulation","query":args.query,"status":"Ambiguous","count":fs.len(),"offset":args.offset,"next_offset":end,"done":end==fs.len(),"formulations":&fs[args.offset..end]})
            }
        },
    };
    out.as_object_mut()
        .unwrap()
        .extend(body.as_object().unwrap().clone());
    out["sources"] = json!(c.edition().sources);
    out["authority"] = json!(
        "Enzyme geometry and formulation documentary evidence are separate. Aliases, isoschizomers and neoschizomers do not transfer supplier conditions or guarantee cleavage."
    );
    Ok(out)
}
pub fn evidence(args: EvidenceQuery) -> Result<Value, ErrorData> {
    let c = builtin()?;
    evidence_for(&c, args)
}
fn evidence_for(c: &Catalog, args: EvidenceQuery) -> Result<Value, ErrorData> {
    check_pin(c, args.catalog_pin.as_ref())?;
    if !bounded(&args.formulation, 128)
        || args.formulation.trim().is_empty()
        || args
            .context
            .as_ref()
            .is_some_and(|s| !bounded(s, 256) || s.trim().is_empty())
    {
        return Err(bad("invalid bounded formulation/context"));
    }
    let condition = match args.condition {
        EvidenceCondition::IncubationTemperature => catalog::Condition::IncubationTemperature,
        EvidenceCondition::BufferActivity {
            buffer,
            temperature_milli_c,
        } => {
            if !bounded(&buffer, 128)
                || buffer.trim().is_empty()
                || !(0..=100000).contains(&temperature_milli_c)
            {
                return Err(bad(
                    "buffer and temperature must describe a bounded valid query",
                ));
            }
            catalog::Condition::BufferActivity {
                buffer,
                temperature_milli_c,
            }
        }
        EvidenceCondition::Methylation { system } => catalog::Condition::Methylation {
            system: match system {
                MethylationSystem::Dam => catalog::MethylationSystem::Dam,
                MethylationSystem::Dcm => catalog::MethylationSystem::Dcm,
                MethylationSystem::CpG => catalog::MethylationSystem::CpG,
            },
        },
    };
    let formulation = match c.lookup_formulation(&args.formulation) {
        Lookup::Unique(f) => f,
        Lookup::Missing => return Err(bad("unknown formulation")),
        Lookup::Ambiguous(_) => {
            return Err(bad("ambiguous formulation; choose an exact catalog ID"));
        }
    };
    let resolution = match c
        .resolve(&formulation.id, &condition, args.context.as_deref())
        .map_err(bad)?
    {
        Resolution::Unknown(evidence) => json!({"status":"Unknown","evidence":evidence}),
        Resolution::Conflict(evidence) => json!({"status":"Conflict","evidence":evidence}),
        Resolution::ContextRequired(contexts) => {
            json!({"status":"ContextRequired","contexts":contexts})
        }
        Resolution::Known { value, evidence } => {
            json!({"status":"Known","value":value,"evidence":evidence})
        }
    };
    Ok(
        json!({"identity":identity(c),"formulation":args.formulation,"formulation_id":formulation.id,"enzyme_id":formulation.enzyme,"formulation_scope":formulation.formulation_scope,"condition":condition,"context":args.context,"resolution":resolution,"sources":c.edition().sources,
        "authority":"Documentary formulation/assay evidence only. Unknown or conflicting information is unresolved; qualifiers and contexts cannot be dropped or generalized into empirical cleavage or buffer compatibility."}),
    )
}

#[cfg(test)]
mod tests {
    use super::{Catalog, EvidenceCondition, EvidenceQuery, builtin, catalog, evidence_for};
    #[test]
    fn serialization_keeps_conflicting_values_and_requires_explicit_assay_context() {
        let c = builtin().unwrap();
        let mut edition = c.edition().clone();
        let f = &mut edition.formulations[0];
        let id = f.id.clone();
        let mut evidence = f
            .evidence
            .iter()
            .find(|e| e.condition == catalog::Condition::IncubationTemperature)
            .unwrap()
            .clone();
        evidence.value = catalog::FactValue::TemperatureMilliC(25000);
        f.evidence.push(evidence.clone());
        let c = Catalog::from_edition(catalog::seal(edition.clone()).unwrap()).unwrap();
        let args = || EvidenceQuery {
            formulation: id.clone(),
            condition: EvidenceCondition::IncubationTemperature,
            context: None,
            catalog_pin: None,
        };
        let conflict = evidence_for(&c, args()).unwrap();
        assert_eq!(conflict["resolution"]["status"], "Conflict");
        assert!(conflict["resolution"].get("value").is_none());
        assert_eq!(
            conflict["resolution"]["evidence"].as_array().unwrap().len(),
            2
        );
        for item in conflict["resolution"]["evidence"].as_array().unwrap() {
            assert!(!item["qualifiers"].as_array().unwrap().is_empty());
            assert!(!item["citations"].as_array().unwrap().is_empty());
        }
        edition.formulations[0].evidence.last_mut().unwrap().context =
            "Independent different assay".into();
        let c = Catalog::from_edition(catalog::seal(edition).unwrap()).unwrap();
        let contexts = evidence_for(&c, args()).unwrap();
        assert_eq!(contexts["resolution"]["status"], "ContextRequired");
        assert_eq!(
            contexts["resolution"]["contexts"].as_array().unwrap().len(),
            2
        );
        let mut explicit = args();
        explicit.context = Some("Independent different assay".into());
        let selected = evidence_for(&c, explicit).unwrap();
        assert_eq!(selected["resolution"]["status"], "Known");
        assert_eq!(
            selected["resolution"]["value"]["temperature_milli_c"],
            25000
        );
    }
}
