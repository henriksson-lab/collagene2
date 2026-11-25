//! Declared oligos from an exact plan, independently of reaction verification.
//! Generated order tables are reference data, never stock, quantities or proof.
use crate::{
    Alphabet, Error, Molecule, Sequence, Topology,
    analysis::Primer,
    chemistry::Terminus,
    formats::ExportReport,
    mutagenesis::PcrOptions,
    oligo_formats::{
        self, ColumnMapping, Delimiter, OligoImportOptions, OligoTable, TableExportMode,
    },
    verification::{self, DesignPlan, Operation, Source},
};
use base64::{Engine, engine::general_purpose::STANDARD};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::Write,
};
const MAX_USES: usize = 2000;
const MAX_REPORT: usize = 32 << 20;
const MAX_ORDERED: usize = 8192;
#[derive(Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct Options {
    pub deduplicate: bool,
}
impl Default for Options {
    fn default() -> Self {
        Self { deduplicate: true }
    }
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OligoUse {
    pub step_index: usize,
    pub kind: String,
    pub primer_name: String,
    pub oligo_id: String,
    pub source: Option<Source>,
    pub primer: Option<Primer>,
    pub mismatch_options: Option<PcrOptions>,
    pub declared_terminus: Option<Terminus>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub plan_revision: String,
    pub target_revision: String,
    pub options: Options,
    pub complete: bool,
    pub unavailable: Vec<String>,
    pub oligos: OligoTable,
    pub uses: Vec<OligoUse>,
    pub warnings: Vec<String>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Exported {
    pub context_revision: String,
    pub report: Report,
    #[serde(with = "export_report")]
    pub export: ExportReport,
}
// Existing text export reports deliberately have no serde implementation. Keep
// the new worker receipt typed without changing other format DTO layouts.
mod export_report {
    use super::ExportReport;
    use serde::{Deserialize, Deserializer, Serialize, Serializer};
    pub fn serialize<S: Serializer>(
        value: &ExportReport,
        serializer: S,
    ) -> Result<S::Ok, S::Error> {
        #[derive(Serialize)]
        struct Wire<'a> {
            text: &'a str,
            warnings: &'a [String],
        }
        Wire {
            text: &value.text,
            warnings: &value.warnings,
        }
        .serialize(serializer)
    }
    pub fn deserialize<'de, D: Deserializer<'de>>(
        deserializer: D,
    ) -> Result<ExportReport, D::Error> {
        #[derive(Deserialize)]
        #[serde(deny_unknown_fields)]
        struct Wire {
            text: String,
            warnings: Vec<String>,
        }
        let value = Wire::deserialize(deserializer)?;
        Ok(ExportReport {
            text: value.text,
            warnings: value.warnings,
        })
    }
}

fn bad(s: &str) -> Error {
    Error::Unsupported(format!("plan oligo output: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|c| c.is_ascii_digit() || (b'a'..=b'f').contains(&c))
}
struct Bounded {
    bytes: Vec<u8>,
    cap: usize,
}
impl Write for Bounded {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        if self
            .bytes
            .len()
            .checked_add(b.len())
            .is_none_or(|n| n > self.cap)
        {
            return Err(std::io::Error::other("bounded plan oligo output"));
        }
        self.bytes.extend_from_slice(b);
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn size(v: &impl Serialize, cap: usize) -> Result<(), Error> {
    struct Count {
        used: usize,
        cap: usize,
    }
    impl Write for Count {
        fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
            self.used = self
                .used
                .checked_add(b.len())
                .filter(|n| *n <= self.cap)
                .ok_or_else(|| std::io::Error::other("bounded plan oligos"))?;
            Ok(b.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    serde_json::to_writer(&mut Count { used: 0, cap }, v).map_err(|_| Error::TooLarge)
}
fn primer_guard(p: &Primer) -> Result<(), Error> {
    if p.name.len() > 1024
        || p.name.chars().any(char::is_control)
        || p.annealing.alphabet() != Alphabet::Dna
        || p.annealing.is_empty()
        || p.annealing
            .len()
            .checked_add(p.tail.len())
            .is_none_or(|n| n > 200)
    {
        return Err(bad(
            "primer name/domain/tail exceeds the typed DNA primer profile",
        ));
    }
    if !p.tail.is_empty() {
        Sequence::dna(&p.tail)?;
    }
    Ok(())
}
fn source_guard(source: &Source, index: usize, ids: &BTreeSet<&str>) -> Result<(), Error> {
    match source {
        Source::Input(id) if id.len() <= 128 && ids.contains(id.as_str()) => Ok(()),
        Source::Step(step) if *step < index => Ok(()),
        _ => Err(bad("unknown, forward or oversized oligo source reference")),
    }
}
/// Borrowed guards precede canonical hashing, source decoding or table allocation.
/// No preparation, annealing, PCR or verification algorithm is invoked.
pub fn preflight(plan: &DesignPlan, target: &Molecule, _options: &Options) -> Result<(), Error> {
    if plan.inputs.len() > 128
        || plan.steps.is_empty()
        || plan.steps.len() > 1000
        || plan.name.len() > 1024
    {
        return Err(Error::TooLarge);
    }
    let mut ids = BTreeSet::new();
    for input in &plan.inputs {
        if input.id.is_empty() || input.id.len() > 128 || !ids.insert(input.id.as_str()) {
            return Err(bad("invalid or duplicate input ID"));
        }
    }
    for (index, step) in plan.steps.iter().enumerate() {
        match step {
            Operation::Pcr {
                template,
                forward,
                reverse,
            }
            | Operation::PcrWithMismatches {
                template,
                forward,
                reverse,
                ..
            } => {
                source_guard(template, index, &ids)?;
                primer_guard(forward)?;
                primer_guard(reverse)?;
            }
            Operation::AnnealOligos { top, bottom, .. } => {
                source_guard(top, index, &ids)?;
                source_guard(bottom, index, &ids)?;
            }
            _ => {}
        }
    }
    crate::assembly_junctions::source_preflight(plan, target)?;
    if target.revision() != plan.target_revision {
        return Err(bad(
            "independent target revision differs from the exact plan binding",
        ));
    }
    for input in &plan.inputs {
        if input.molecule.revision() != input.revision {
            return Err(bad("original input revision changed"));
        }
    }
    Ok(())
}
fn warning(warnings: &mut Vec<String>, message: String) {
    if warnings.len() < 128 {
        let mut text = message.chars().take(2048).collect::<String>();
        if text.len() < message.len() {
            text.push('…');
        }
        warnings.push(text);
    } else if warnings.len() == 128 {
        warnings.push("Further order-output warnings omitted; inspect the full typed plan and original sources.".into());
    }
}
struct Entry {
    name: String,
    ordered: String,
    domains: Option<(String, String)>,
    modification: Option<String>,
    names: Vec<String>,
}
struct Builder {
    entries: Vec<Entry>,
    uses: Vec<(usize, OligoUse)>,
    keys: BTreeMap<(String, Option<String>), usize>,
    warnings: Vec<String>,
    unavailable: Vec<String>,
    dedup: bool,
}
impl Builder {
    fn push(
        &mut self,
        name: String,
        ordered: String,
        domains: Option<(String, String)>,
        modification: Option<String>,
        usage: OligoUse,
    ) -> Result<(), Error> {
        if self.uses.len() >= MAX_USES || ordered.is_empty() || ordered.len() > MAX_ORDERED {
            return Err(Error::TooLarge);
        }
        let key = (ordered.clone(), modification.clone());
        let existing = if self.dedup {
            self.keys.get(&key).copied()
        } else {
            None
        };
        let index = if let Some(index) = existing {
            let entry = &mut self.entries[index];
            if entry.domains != domains {
                entry.domains = None;
                warning(
                    &mut self.warnings,
                    format!(
                        "Order row {} has multiple domain/tail partitions or an unassigned anneal domain. Domain columns are blank; every original primer partition remains in uses.",
                        index + 1
                    ),
                );
            }
            if !entry.names.contains(&name) {
                warning(
                    &mut self.warnings,
                    format!(
                        "Order row {} shares ordered DNA/modifications across distinct names. Original names remain separate in uses; the order-row name is only the first label.",
                        index + 1
                    ),
                );
                entry.names.push(name);
            }
            index
        } else {
            let index = self.entries.len();
            self.keys.entry(key).or_insert(index);
            self.entries.push(Entry {
                name: name.clone(),
                ordered,
                domains,
                modification,
                names: vec![name],
            });
            index
        };
        self.uses.push((index, usage));
        Ok(())
    }
    fn pcr(
        &mut self,
        index: usize,
        kind: &str,
        source: &Source,
        p: &Primer,
        policy: Option<&PcrOptions>,
    ) -> Result<(), Error> {
        let normalized = Primer::new(&p.name, p.annealing.as_str(), &p.tail)?;
        let name = if p.name.is_empty() {
            warning(
                &mut self.warnings,
                format!(
                    "Step {} {kind} primer has no name. A generated order-row label is used; the original empty name remains in uses.",
                    index + 1
                ),
            );
            format!("unnamed-step-{}-{kind}", index + 1)
        } else {
            p.name.clone()
        };
        self.push(
            name,
            normalized.ordered_sequence(),
            Some((normalized.annealing.as_str().into(), normalized.tail)),
            None,
            OligoUse {
                step_index: index,
                kind: kind.into(),
                primer_name: p.name.clone(),
                oligo_id: String::new(),
                source: Some(source.clone()),
                primer: Some(p.clone()),
                mismatch_options: policy.cloned(),
                declared_terminus: None,
            },
        )
    }
    fn anneal(
        &mut self,
        plan: &DesignPlan,
        index: usize,
        kind: &str,
        source: &Source,
        terminus: Terminus,
    ) -> Result<(), Error> {
        let Source::Input(id) = source else {
            warning(
                &mut self.warnings,
                format!(
                    "Step {} {kind}: ordered oligo unavailable for a derived Step source. Reconstruct its originating graph; no descendant sequence was guessed as an order.",
                    index + 1
                ),
            );
            return Ok(());
        };
        let input = plan
            .inputs
            .iter()
            .find(|input| &input.id == id)
            .ok_or_else(|| bad("missing original oligo"))?;
        let molecule = &input.molecule;
        if molecule.sequence.alphabet() != Alphabet::Dna
            || molecule.topology != Topology::Linear
            || molecule.sequence.len() > MAX_ORDERED
        {
            warning(
                &mut self.unavailable,
                format!(
                    "Step {} {kind}: original source is not a bounded linear DNA oligo; no order row created.",
                    index + 1
                ),
            );
            return Ok(());
        }
        let modification = match source_modification(molecule) {
            Ok(value) => value,
            Err(error) => {
                warning(
                    &mut self.unavailable,
                    format!(
                        "Step {} {kind}: order unavailable: {error}. Retain original source fields; modification-bearing evidence was not flattened.",
                        index + 1
                    ),
                );
                return Ok(());
            }
        };
        warning(
            &mut self.warnings,
            format!(
                "Step {} {kind}: declared termini are state assumptions only, not supplier modification requests. No annealing-domain partition or order quantity was inferred.",
                index + 1
            ),
        );
        let name = if molecule.name.is_empty() {
            format!("unnamed-step-{}-{kind}", index + 1)
        } else {
            molecule.name.clone()
        };
        if name.len() > 1024 || name.chars().any(char::is_control) {
            return Err(bad("original oligo name is not a supported order label"));
        }
        self.push(
            name,
            molecule.sequence.as_str().into(),
            None,
            modification,
            OligoUse {
                step_index: index,
                kind: kind.into(),
                primer_name: molecule.name.clone(),
                oligo_id: String::new(),
                source: Some(source.clone()),
                primer: None,
                mismatch_options: None,
                declared_terminus: Some(terminus),
            },
        )
    }
}
/// Only existing parser-backed original table evidence may supply an explicit
/// modification string. Metadata is never interpreted as current chemistry.
fn source_modification(m: &Molecule) -> Result<Option<String>, Error> {
    if !m.metadata.keys().any(|key| key.starts_with("oligo:")) {
        return Ok(None);
    }
    let raw = m
        .metadata
        .get("oligo:original_record_base64")
        .ok_or_else(|| bad("incomplete original oligo record"))?;
    if raw.len() > 128 * 1024 {
        return Err(Error::TooLarge);
    }
    let bytes = STANDARD
        .decode(raw)
        .map_err(|_| bad("invalid original oligo record encoding"))?;
    let expected = m
        .metadata
        .get("oligo:original_record_blake3")
        .ok_or_else(|| bad("missing original oligo record hash"))?;
    if blake3::hash(&bytes).to_hex().as_str() != expected {
        return Err(bad("original oligo record hash differs"));
    }
    let mapping = m
        .metadata
        .get("oligo:source_mapping")
        .filter(|s| s.len() <= 8192)
        .ok_or_else(|| bad("missing bounded original column mapping"))?;
    let mut options: OligoImportOptions =
        serde_json::from_str(mapping).map_err(|_| bad("invalid original column mapping"))?;
    options.has_header = false;
    let text =
        std::str::from_utf8(&bytes).map_err(|_| bad("original oligo record is not UTF-8"))?;
    let table = oligo_formats::import_oligos(text, &options)?;
    if table.records.len() != 1 {
        return Err(bad(
            "original oligo evidence is not one complete source row",
        ));
    }
    let record = &table.records[0];
    let ordered = record
        .raw_ordered_sequence
        .as_deref()
        .map(Sequence::dna)
        .transpose()?
        .or_else(|| {
            record
                .primer
                .as_ref()
                .and_then(|p| Sequence::dna(&p.ordered_sequence()).ok())
        })
        .ok_or_else(|| bad("no faithful literal ordered DNA representation"))?;
    if ordered != m.sequence {
        return Err(bad("original oligo order differs from current source DNA"));
    }
    if let Some(fields) = m.metadata.get("oligo:source_fields") {
        if fields.len() > 128 * 1024 {
            return Err(Error::TooLarge);
        }
        let fields: Vec<String> =
            serde_json::from_str(fields).map_err(|_| bad("invalid original fields"))?;
        if fields != record.original_fields {
            return Err(bad("original oligo field evidence differs"));
        }
    }
    if record.modification.is_some() {
        return Err(bad(
            "explicit original modification notation is opaque to the current parser; preserve the original record and choose a faithful supplier representation explicitly",
        ));
    }
    Ok(None)
}
pub fn build(plan: &DesignPlan, target: &Molecule, options: &Options) -> Result<Report, Error> {
    preflight(plan, target, options)?;
    let mut builder=Builder{entries:vec![],uses:vec![],keys:BTreeMap::new(),warnings:vec!["Declared oligos from the exact typed plan only. Order extraction performs no PCR or target verification and does not establish experimental success.".into(),"No stock, supplier notation, order quantities or downstream treatment-derived modifications are inferred. PCR primers have no explicit modification field in this model.".into()],unavailable:vec![],dedup:options.deduplicate};
    for (index, step) in plan.steps.iter().enumerate() {
        match step {
            Operation::Pcr {
                template,
                forward,
                reverse,
            } => {
                builder.pcr(index, "pcr-forward", template, forward, None)?;
                builder.pcr(index, "pcr-reverse", template, reverse, None)?;
            }
            Operation::PcrWithMismatches {
                template,
                forward,
                reverse,
                options,
            } => {
                builder.pcr(
                    index,
                    "mismatch-pcr-forward",
                    template,
                    forward,
                    Some(options),
                )?;
                builder.pcr(
                    index,
                    "mismatch-pcr-reverse",
                    template,
                    reverse,
                    Some(options),
                )?;
            }
            Operation::AnnealOligos {
                top,
                bottom,
                options,
                ..
            } => {
                builder.anneal(plan, index, "anneal-top", top, options.top_chemistry)?;
                builder.anneal(
                    plan,
                    index,
                    "anneal-bottom",
                    bottom,
                    options.bottom_chemistry,
                )?;
            }
            Operation::MultiSiteRecombineReference { .. } => warning(
                &mut builder.warnings,
                format!(
                    "Step {} MultiSite reference exchange supplies no PCR oligos. No attachment-site primer or upstream entry-clone preparation was inferred; only separately declared typed PCR steps contribute primers.",
                    index + 1
                ),
            ),
            Operation::Topo { .. }
            | Operation::CheckedTopo { .. }
            | Operation::ProcessedTopo { .. }
            | Operation::TopoActions { .. } => warning(
                &mut builder.unavailable,
                format!(
                    "Step {} embedded TOPO preparation oligos are outside this typed-PCR extraction profile; inspect the retained physical recipe separately.",
                    index + 1
                ),
            ),
            _ => {}
        }
    }
    if builder.entries.is_empty() {
        warning(&mut builder.warnings,"No orderable typed PCR primers or original anneal oligos were found. No synthesis order was inferred from other preparation steps.".into());
    }
    let mut writer = csv::WriterBuilder::new()
        .terminator(csv::Terminator::CRLF)
        .from_writer(Bounded {
            bytes: vec![],
            cap: oligo_formats::MAX_TABLE_BYTES,
        });
    writer
        .write_record([
            "Name",
            "Ordered sequence",
            "Annealing domain",
            "5prime tail",
            "Modifications",
            "Names preview",
            "Distinct labels (not quantity)",
        ])
        .map_err(|_| Error::TooLarge)?;
    for entry in &builder.entries {
        let (domain, tail) = entry
            .domains
            .as_ref()
            .map(|(d, t)| (d.as_str(), t.as_str()))
            .unwrap_or(("", ""));
        let names = entry
            .names
            .iter()
            .take(16)
            .map(|name| name.chars().take(64).collect::<String>())
            .collect::<Vec<_>>()
            .join(" | ");
        writer
            .write_record([
                entry.name.as_str(),
                entry.ordered.as_str(),
                domain,
                tail,
                entry.modification.as_deref().unwrap_or(""),
                &names,
                &entry.names.len().to_string(),
            ])
            .map_err(|_| Error::TooLarge)?;
    }
    let bytes = writer.into_inner().map_err(|_| Error::TooLarge)?.bytes;
    let text = String::from_utf8(bytes).map_err(|_| bad("generated order table is not UTF-8"))?;
    let oligos = oligo_formats::import_oligos(
        &text,
        &OligoImportOptions {
            delimiter: Delimiter::Comma,
            has_header: true,
            columns: ColumnMapping {
                name: 0,
                ordered_sequence: Some(1),
                annealing_sequence: Some(2),
                tail: Some(3),
                modification: Some(4),
                amount: None,
                unit: None,
            },
        },
    )?;
    let uses = builder
        .uses
        .into_iter()
        .map(|(row, mut usage)| {
            usage.oligo_id = oligos.records[row].row_id.clone();
            usage
        })
        .collect();
    let report = Report {
        plan_revision: verification::plan_revision(plan),
        target_revision: target.revision(),
        options: options.clone(),
        complete: builder.unavailable.is_empty(),
        unavailable: builder.unavailable,
        oligos,
        uses,
        warnings: builder.warnings,
    };
    report.validate()?;
    Ok(report)
}
impl Report {
    pub fn validate(&self) -> Result<(), Error> {
        if !hex(&self.plan_revision)
            || !hex(&self.target_revision)
            || self.uses.len() > MAX_USES
            || self.oligos.records.len() > MAX_USES
            || self.oligos.original_base64.len() > oligo_formats::MAX_TABLE_BYTES.div_ceil(3) * 4
            || self.warnings.len() > 129
            || self.unavailable.len() > 129
            || self.complete != self.unavailable.is_empty()
        {
            return Err(Error::TooLarge);
        }
        for usage in &self.uses {
            if usage.kind.len() > 64
                || usage.primer_name.len() > 1024
                || usage.oligo_id.len() > 96
                || usage
                    .oligo_id
                    .split_once(':')
                    .is_none_or(|(hash, index)| !hex(hash) || index.parse::<usize>().is_err())
                || usage.step_index >= 1000
            {
                return Err(bad("invalid use identity"));
            }
            if let Some(primer) = &usage.primer {
                primer_guard(primer)?;
            }
            if let Some(Source::Input(id)) = &usage.source {
                if id.len() > 128 {
                    return Err(Error::TooLarge);
                }
            }
        }
        if self.oligos.headers.len() > 64
            || self.oligos.headers.iter().any(|s| s.len() > 8192)
            || !hex(&self.oligos.source_blake3)
            || self.oligos.warnings.len() > 129
        {
            return Err(Error::TooLarge);
        }
        for warning in self
            .warnings
            .iter()
            .chain(&self.unavailable)
            .chain(&self.oligos.warnings)
        {
            if warning.len() > 8196 {
                return Err(Error::TooLarge);
            }
        }
        for row in &self.oligos.records {
            if row.name.len() > 1024
                || row.row_id.len() > 96
                || !hex(&row.original_record_blake3)
                || row.original_fields.len() > 64
                || row.original_fields.iter().any(|field| field.len() > 8192)
                || row.original_record_base64.len() > 128 * 1024
                || row.quantity.is_some()
                || row.diagnostics.len() > 32
                || row.diagnostics.iter().any(|s| s.len() > 8196)
                || row
                    .raw_ordered_sequence
                    .as_ref()
                    .is_some_and(|s| s.len() > MAX_ORDERED)
                || row
                    .raw_annealing_sequence
                    .as_ref()
                    .is_some_and(|s| s.len() > MAX_ORDERED)
                || row.raw_tail.as_ref().is_some_and(|s| s.len() > MAX_ORDERED)
                || row.modification.as_ref().is_some_and(|s| s.len() > 8192)
            {
                return Err(Error::TooLarge);
            }
            if let Some(primer) = &row.primer {
                primer_guard(primer)?;
            }
        }
        size(self, MAX_REPORT)?;
        self.oligos.validate()?;
        let ids = self
            .oligos
            .records
            .iter()
            .map(|r| r.row_id.as_str())
            .collect::<BTreeSet<_>>();
        if self
            .uses
            .iter()
            .any(|usage| !ids.contains(usage.oligo_id.as_str()))
        {
            return Err(bad("order usage does not name a retained row"));
        }
        Ok(())
    }
}
pub fn export_revision(
    plan: &DesignPlan,
    target: &Molecule,
    options: &Options,
    delimiter: Delimiter,
) -> Result<String, Error> {
    preflight(plan, target, options)?;
    crate::assembly_recipe::hash(&(
        "plan-oligo-export-1",
        verification::plan_revision(plan),
        target.revision(),
        options,
        delimiter,
    ))
}
pub fn export(
    plan: &DesignPlan,
    target: &Molecule,
    options: &Options,
    delimiter: Delimiter,
) -> Result<Exported, Error> {
    let context_revision = export_revision(plan, target, options, delimiter)?;
    let report = build(plan, target, options)?;
    let export =
        oligo_formats::export_oligos(&report.oligos, delimiter, TableExportMode::SpreadsheetSafe)?;
    let value = Exported {
        context_revision,
        report,
        export,
    };
    value.validate_binding(&value.context_revision)?;
    Ok(value)
}
impl Exported {
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if !hex(expected)
            || self.context_revision != expected
            || self.export.text.len() > oligo_formats::MAX_TABLE_BYTES
            || self.export.warnings.len() > 129
            || self.export.warnings.iter().any(|s| s.len() > 8196)
        {
            return Err(bad("export context or bounded text differs"));
        }
        self.report.validate()?;
        size(self, MAX_REPORT)
    }
}
