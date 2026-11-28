//! Published, symmetric ligation-frequency weights from Pryor et al. (2020).
//! Original third-party data is CC-BY-4.0; this independent importer/scorer is MIT.
//! These are published pair-frequency entries, not independent raw directional
//! sequencing observations. No optimizer implementation or Excel runtime is used.
use crate::cloning::{
    assembly::recipe as assembly_recipe, golden_gate::GoldenGateEnzyme,
    ligation::fidelity as numeric,
};
use collagene_model::Error;
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet},
    io::{Cursor, Read},
    sync::{Arc, OnceLock},
};
pub const MODEL: &str = "published-symmetric-pair-frequencies-1";
pub const VERSION: &str = "pryor2020-plos-s1-s5-1";
pub const MAX_SOURCE_BYTES: usize = 1 << 20;
const MAX_INFLATED: usize = 8 << 20;
const MAX_XML: usize = 4 << 20;
const NS: &str = "http://schemas.openxmlformats.org/spreadsheetml/2006/main";
macro_rules! dto {($($t:tt)*)=>{#[derive(Clone,Debug,PartialEq,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct ProfileRef{pub profile_id:String,pub profile_revision:String}}
dto! {pub struct Metadata{
    pub profile:ProfileRef,pub version:String,pub model:String,pub name:String,
    /// Geometry family only; exact enzyme formulation is a separate declaration.
    pub geometry_family:GoldenGateEnzyme,pub enzyme_formulation:String,pub overhang_length:usize,
    pub conditions:numeric::Conditions,pub provenance:numeric::Provenance,pub source_sha256:String,pub cells:usize,pub data_revision:String,
}}
dto! {pub struct Summary{
    pub profile:ProfileRef,pub version:String,pub model:String,pub name:String,
    pub geometry_family:GoldenGateEnzyme,pub enzyme_formulation:String,pub overhang_length:usize,
    pub conditions:numeric::Conditions,pub provenance:numeric::Provenance,pub source_sha256:String,pub cells:usize,
}}
dto! {pub struct Orientation{
    pub end_5to3:String,pub correct_partner_5to3:String,pub correct_weight:Option<u64>,
    #[serde(with="numeric::optional_decimal_u128")]
    #[cfg_attr(feature="json-schema",schemars(with="Option<String>"))]
    pub total_weight:Option<u128>,pub missing_partners:Vec<String>,
}}
dto! {pub enum Estimate{
    Known{
        #[serde(with="numeric::decimal_u128")]
        #[cfg_attr(feature="json-schema",schemars(with="String"))]
        correct_weight:u128,
        #[serde(with="numeric::decimal_u128")]
        #[cfg_attr(feature="json-schema",schemars(with="String"))]
        total_weight:u128,fidelity:f64,
    },Unknown{reason:numeric::UnknownReason},
}}
dto! {pub struct PairScore{pub canonical_end:String,pub reverse_complement:String,pub occurrences:Vec<usize>,pub orientations:Vec<Orientation>,pub estimate:Estimate}}
dto! {pub struct Report{pub metadata:Metadata,pub context_revision:String,pub requested_ends:Vec<String>,pub physical_competitors:Vec<String>,pub physical_ends:Vec<String>,pub pairs:Vec<PairScore>,pub overall_fidelity:Option<f64>,pub log_overall_fidelity:Option<f64>,pub warnings:Vec<String>}}
fn bad(s: &str) -> Error {
    Error::Format(format!("published ligation profile: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn index(s: &str, k: usize) -> Option<usize> {
    if s.len() != k {
        return None;
    }
    s.bytes().try_fold(0usize, |n, b| {
        Some(
            n * 4
                + match b {
                    b'A' => 0,
                    b'C' => 1,
                    b'G' => 2,
                    b'T' => 3,
                    _ => return None,
                },
        )
    })
}
fn xml(bytes: &[u8]) -> Result<roxmltree::Document<'_>, Error> {
    if bytes.len() > MAX_XML {
        return Err(Error::TooLarge);
    }
    let text = std::str::from_utf8(bytes).map_err(|_| bad("XML must be UTF8"))?;
    if text.contains("<!DOCTYPE") || text.contains("<!ENTITY") {
        return Err(bad("DTD/entities are forbidden"));
    }
    roxmltree::Document::parse_with_options(
        text,
        roxmltree::ParsingOptions {
            allow_dtd: false,
            nodes_limit: 300_000,
            ..Default::default()
        },
    )
    .map_err(|_| bad("invalid/bounded XML"))
}
fn coordinate(text: &str, max: usize) -> Result<(usize, usize), Error> {
    if text.len() > 8 {
        return Err(bad("cell coordinate bound"));
    }
    let mut column = 0usize;
    let mut at = 0;
    for b in text.bytes().take_while(|b| b.is_ascii_uppercase()) {
        column = column
            .checked_mul(26)
            .and_then(|n| n.checked_add((b - b'A' + 1) as usize))
            .ok_or(Error::TooLarge)?;
        at += 1;
    }
    let row = text[at..]
        .parse::<usize>()
        .map_err(|_| bad("invalid literal cell coordinate"))?;
    if column == 0 || column > max || row == 0 || row > max || at == 0 {
        return Err(bad("cell outside complete matrix"));
    }
    Ok((row - 1, column - 1))
}
/// Bounded independently written XLSX matrix reader. The result has no profile
/// identity/license authority; only load() installs an exact embedded source.
pub struct Matrix {
    overhang_length: usize,
    weights: Vec<u64>,
    source_blake3: String,
}
impl Matrix {
    pub fn overhang_length(&self) -> usize {
        self.overhang_length
    }
    pub fn cell_count(&self) -> usize {
        self.weights.len()
    }
    pub fn source_blake3(&self) -> &str {
        &self.source_blake3
    }
    pub fn weight(&self, end: &str, partner: &str) -> Option<u64> {
        let n = 1usize << (2 * self.overhang_length);
        Some(
            *self.weights.get(
                index(end, self.overhang_length)? * n + index(partner, self.overhang_length)?,
            )?,
        )
    }
}
pub fn parse_published_matrix(bytes: &[u8], k: usize) -> Result<Matrix, Error> {
    if bytes.len() > MAX_SOURCE_BYTES || !matches!(k, 3 | 4) {
        return Err(Error::TooLarge);
    }
    let mut zip = zip::ZipArchive::new(Cursor::new(bytes)).map_err(|_| bad("invalid XLSX ZIP"))?;
    const PARTS: [&str; 10] = [
        "[Content_Types].xml",
        "_rels/.rels",
        "xl/workbook.xml",
        "xl/_rels/workbook.xml.rels",
        "xl/worksheets/sheet1.xml",
        "xl/theme/theme1.xml",
        "xl/styles.xml",
        "xl/sharedStrings.xml",
        "docProps/core.xml",
        "docProps/app.xml",
    ];
    if zip.len() != PARTS.len() {
        return Err(bad("unexpected/missing workbook parts"));
    }
    let mut parts = BTreeMap::<String, Vec<u8>>::new();
    let mut inflated = 0usize;
    for i in 0..zip.len() {
        let mut entry = zip
            .by_index(i)
            .map_err(|_| bad("unreadable workbook entry"))?;
        let name = entry.name().to_owned();
        if !PARTS.contains(&name.as_str())
            || parts.contains_key(&name)
            || entry.is_dir()
            || entry.enclosed_name().is_none()
            || entry.size() > MAX_XML as u64
        {
            return Err(bad("external/duplicate/oversized workbook part"));
        }
        inflated = inflated
            .checked_add(usize::try_from(entry.size()).map_err(|_| Error::TooLarge)?)
            .filter(|n| *n <= MAX_INFLATED)
            .ok_or(Error::TooLarge)?;
        let mut raw = Vec::with_capacity(entry.size() as usize);
        entry
            .by_ref()
            .take((MAX_XML + 1) as u64)
            .read_to_end(&mut raw)
            .map_err(|_| bad("ZIP integrity/read failure"))?;
        if raw.len() > MAX_XML || raw.len() as u64 != entry.size() {
            return Err(bad("truncated/oversized inflated XML"));
        }
        let doc = xml(&raw)?;
        if doc
            .descendants()
            .filter(|n| n.is_element())
            .any(|n| n.attribute("TargetMode") == Some("External") || n.tag_name().name() == "f")
        {
            return Err(bad("external relationship or formula is forbidden"));
        }
        if name.ends_with(".rels")
            && doc
                .descendants()
                .filter(|n| n.is_element() && n.tag_name().name() == "Relationship")
                .any(|n| {
                    n.attribute("Target")
                        .is_none_or(|t| t.contains("://") || t.starts_with('/') || t.contains(".."))
                })
        {
            return Err(bad("external/relative-escape relationship"));
        }
        if name.ends_with(".rels") {
            let mut ids = BTreeSet::new();
            for relation in doc
                .descendants()
                .filter(|n| n.is_element() && n.tag_name().name() == "Relationship")
            {
                let id = relation
                    .attribute("Id")
                    .ok_or_else(|| bad("relationship ID missing"))?;
                if id.len() > 64 || !ids.insert(id) {
                    return Err(bad("duplicate/oversized relationship ID"));
                }
            }
        }
        parts.insert(name, raw);
    }
    let book = xml(&parts["xl/workbook.xml"])?;
    let sheets = book
        .descendants()
        .filter(|n| n.has_tag_name((NS, "sheet")))
        .collect::<Vec<_>>();
    if sheets.len() != 1
        || sheets[0].attribute((
            "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
            "id",
        )) != Some("rId1")
    {
        return Err(bad("one fixed worksheet required"));
    }
    let relationships = xml(&parts["xl/_rels/workbook.xml.rels"])?;
    if !relationships.descendants().any(|n| {
        n.attribute("Id") == Some("rId1") && n.attribute("Target") == Some("worksheets/sheet1.xml")
    }) {
        return Err(bad("worksheet relationship changed"));
    }
    let strings_doc = xml(&parts["xl/sharedStrings.xml"])?;
    let universe = 1usize << (2 * k);
    let mut strings = vec![];
    for si in strings_doc
        .root_element()
        .children()
        .filter(|n| n.has_tag_name((NS, "si")))
    {
        if strings.len() > universe * 2 {
            return Err(Error::TooLarge);
        }
        let value = si
            .descendants()
            .filter(|n| n.has_tag_name((NS, "t")))
            .filter_map(|n| n.text())
            .collect::<String>();
        if value.len() > 32 {
            return Err(bad("shared string bound"));
        }
        strings.push(value);
    }
    let sheet = xml(&parts["xl/worksheets/sheet1.xml"])?;
    if !sheet.root_element().has_tag_name((NS, "worksheet")) {
        return Err(bad("invalid worksheet namespace/root"));
    }
    let sheetdata = sheet
        .root_element()
        .children()
        .filter(|n| n.has_tag_name((NS, "sheetData")))
        .collect::<Vec<_>>();
    if sheetdata.len() != 1 {
        return Err(bad("exactly one sheetData required"));
    }
    let width = universe + 1;
    let mut cells = vec![None::<String>; width * width];
    for row in sheetdata[0]
        .children()
        .filter(|n| n.has_tag_name((NS, "row")))
    {
        for cell in row.children().filter(|n| n.has_tag_name((NS, "c"))) {
            let (r, c) = coordinate(
                cell.attribute("r")
                    .ok_or_else(|| bad("missing cell coordinate"))?,
                width,
            )?;
            if cells[r * width + c].is_some() {
                return Err(bad("duplicate cell"));
            }
            let values = cell
                .children()
                .filter(|n| n.has_tag_name((NS, "v")))
                .collect::<Vec<_>>();
            if values.len() != 1 {
                return Err(bad(
                    "literal single value required; formulas/missing values unsupported",
                ));
            }
            let raw = values[0].text().ok_or_else(|| bad("empty cell value"))?;
            if raw.len() > 32 {
                return Err(bad("bounded cell value required"));
            }
            let value = match cell.attribute("t") {
                Some("s") => strings
                    .get(
                        raw.parse::<usize>()
                            .map_err(|_| bad("invalid shared index"))?,
                    )
                    .cloned()
                    .ok_or_else(|| bad("shared string index missing"))?,
                None | Some("n") => raw.into(),
                _ => return Err(bad("only literal shared strings/integer cells supported")),
            };
            if r == 0 || c == 0 {
                if cell.attribute("t") != Some("s") {
                    return Err(bad("matrix headers must be literal end labels"));
                }
            } else if cell.attribute("t") == Some("s") {
                return Err(bad("frequency values must be numeric, not strings"));
            }
            cells[r * width + c] = Some(value);
        }
    }
    if cells.iter().any(Option::is_none) || cells[0].as_deref() != Some("Overhang") {
        return Err(bad("incomplete matrix or changed header"));
    }
    let mut rows = vec![];
    let mut columns = vec![];
    let mut seen_rows = BTreeSet::new();
    let mut seen_columns = BTreeSet::new();
    for i in 1..width {
        let r = cells[i * width].as_deref().ok_or(Error::TooLarge)?;
        let c = cells[i].as_deref().ok_or(Error::TooLarge)?;
        let r = index(r, k).ok_or_else(|| bad("invalid row end label"))?;
        let c = index(c, k).ok_or_else(|| bad("invalid column end label"))?;
        if !seen_rows.insert(r) || !seen_columns.insert(c) {
            return Err(bad("duplicate physical end label"));
        }
        rows.push(r);
        columns.push(c);
    }
    let mut weights = vec![0u64; universe * universe];
    for (r, ri) in rows.iter().enumerate() {
        for (c, ci) in columns.iter().enumerate() {
            let raw = cells[(r + 1) * width + c + 1]
                .as_deref()
                .ok_or(Error::TooLarge)?;
            if raw.is_empty()
                || !raw.bytes().all(|b| b.is_ascii_digit())
                || (raw.len() > 1 && raw.starts_with('0'))
            {
                return Err(bad("frequency must be an exact nonnegative integer"));
            }
            weights[ri * universe + ci] = raw
                .parse::<u64>()
                .map_err(|_| bad("frequency weight overflow"))?;
        }
    }
    for r in 0..universe {
        for c in 0..universe {
            if weights[r * universe + c] != weights[c * universe + r] {
                return Err(bad("published symmetric pair-frequency invariant changed"));
            }
        }
    }
    Ok(Matrix {
        overhang_length: k,
        weights,
        source_blake3: blake3::hash(bytes).to_hex().to_string(),
    })
}
struct SourceDef {
    id: &'static str,
    enzyme: &'static str,
    k: usize,
    family: GoldenGateEnzyme,
    sha256: &'static str,
    table: usize,
    raw: &'static [u8],
}
static SOURCES: [SourceDef; 5] = [
    SourceDef {
        id: "pryor2020-bsai-hfv2",
        enzyme: "BsaI-HFv2",
        k: 4,
        family: GoldenGateEnzyme::BsaI,
        sha256: "320dd058f3ca6372768c7ddfe9cda2a2ea2a076c4f86dd08587d8f295aae1d15",
        table: 1,
        raw: include_bytes!("../../../../../data/pryor2020/pone.0238592.s001.xlsx"),
    },
    SourceDef {
        id: "pryor2020-bsmbi-v2",
        enzyme: "BsmBI-v2",
        k: 4,
        family: GoldenGateEnzyme::BsmBI,
        sha256: "7c444e99e5e4d245461a892e8543a35c84987f93abcee7098de9d9d817237e94",
        table: 2,
        raw: include_bytes!("../../../../../data/pryor2020/pone.0238592.s002.xlsx"),
    },
    SourceDef {
        id: "pryor2020-esp3i",
        enzyme: "Esp3I",
        k: 4,
        family: GoldenGateEnzyme::BsmBI,
        sha256: "1557e62cfa4f89cd021420b75e145dae2f453ecf26dcca340a8c29008736ea66",
        table: 3,
        raw: include_bytes!("../../../../../data/pryor2020/pone.0238592.s003.xlsx"),
    },
    SourceDef {
        id: "pryor2020-bbsi-hf",
        enzyme: "BbsI-HF",
        k: 4,
        family: GoldenGateEnzyme::BbsI,
        sha256: "56143bb445e6d84bba646429402e0948793395269cc5d8ba0e80962c0cac6492",
        table: 4,
        raw: include_bytes!("../../../../../data/pryor2020/pone.0238592.s004.xlsx"),
    },
    SourceDef {
        id: "pryor2020-sapi",
        enzyme: "SapI",
        k: 3,
        family: GoldenGateEnzyme::SapI,
        sha256: "dfa2df906bd306e3971652c4b95d004818fdb8dd9e49f4522691959d28f537cb",
        table: 5,
        raw: include_bytes!("../../../../../data/pryor2020/pone.0238592.s005.xlsx"),
    },
];
fn summary(def: &SourceDef) -> Result<Summary, Error> {
    let high = if def.table == 2 { 42 } else { 37 };
    let formulation = if def.table <= 2 {
        format!(
            "2 µL corresponding {} Golden Gate Enzyme Mix per20 µL assay; component enzyme units not separately specified",
            def.enzyme
        )
    } else {
        format!("500 U T4 DNA ligase +15 U {} per20 µL assay", def.enzyme)
    };
    let conditions=numeric::Conditions{ligase:"T4 DNA ligase".into(),buffer:if def.table==4{"1×CutSmart buffer supplemented with10 mM DTT and1 mM ATP"}else{"1×T4 DNA ligase buffer"}.into(),temperature_protocol:format!("30 cycles:5 min at{high}°C then5 min at16°C; final5 min at60°C"),incubation:format!("300 min cycling plus5 min final heat step; {formulation}"),restriction_enzymes:vec![def.enzyme.into()],substrate_context:"100 nM randomized-overhang hairpin assay substrate,20 µL reaction; sequencing of ligated products; at least two experiments on different days. Table does not identify separate replicate values or their pooling/normalization.".into()};
    let provenance=numeric::Provenance{source_url:format!("https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0238592.s{:03}&type=supplementary",def.table),citation:format!("Pryor JM, Potapov V, Kucera RB, Bilotti K, Cantor EJ, Lohman GJS (2020). Enabling one-pot Golden Gate assemblies of unprecedented complexity using data-optimized assembly design. PLOS ONE15(9):e0238592. doi:10.1371/journal.pone.0238592; S{} Table.",def.table),license:"CC-BY-4.0; https://creativecommons.org/licenses/by/4.0/".into(),attribution:"©2020 Pryor et al. Published by PLOS. Data separately licensed from MIT importer/scorer; no NEB/PLOS endorsement implied.".into(),source_blake3:blake3::hash(def.raw).to_hex().to_string(),normalization:format!("{VERSION}: literal integer published ligation-frequency entries, mapped by actual row/column labels. Published table is symmetric; original replicate/sequence-read preprocessing unspecified. No transposition fill, symmetry averaging, rescaling, row normalization or fabricated independent directional observations.")};
    let profile_revision = assembly_recipe::hash(&(
        MODEL,
        VERSION,
        def.id,
        def.enzyme,
        def.k,
        def.family,
        &conditions,
        &provenance,
        def.sha256,
    ))?;
    Ok(Summary {
        profile: ProfileRef {
            profile_id: def.id.into(),
            profile_revision,
        },
        version: VERSION.into(),
        model: MODEL.into(),
        name: format!("Pryor2020 · {} · T4 ligase · {high}/16°C", def.enzyme),
        geometry_family: def.family,
        enzyme_formulation: def.enzyme.into(),
        overhang_length: def.k,
        conditions,
        provenance,
        source_sha256: def.sha256.into(),
        cells: 1usize << (4 * def.k),
    })
}
static SUMMARIES: OnceLock<Result<Vec<Summary>, String>> = OnceLock::new();
static PROFILES: [OnceLock<Result<Arc<Profile>, String>>; 5] = [const { OnceLock::new() }; 5];
/// Fixed metadata and hashes only: no XLSX parsing or scientific enumeration.
pub fn list() -> Result<Vec<Summary>, Error> {
    SUMMARIES
        .get_or_init(|| {
            SOURCES
                .iter()
                .map(summary)
                .collect::<Result<Vec<_>, _>>()
                .map_err(|e| e.to_string())
        })
        .as_ref()
        .cloned()
        .map_err(|e| bad(e))
}
pub fn summaries() -> Result<Vec<Summary>, Error> {
    list()
}
pub fn catalog_revision() -> Result<String, Error> {
    assembly_recipe::hash(&(VERSION, list()?))
}
fn resolve(reference: &ProfileRef) -> Result<(usize, Summary), Error> {
    if reference.profile_id.len() > 64 || !hex(&reference.profile_revision) {
        return Err(bad("bounded exact profile reference required"));
    }
    let all = list()?;
    all.into_iter()
        .enumerate()
        .find(|(_, s)| s.profile == *reference)
        .ok_or_else(|| bad("published profile identity/version/source changed"))
}
pub struct Profile {
    metadata: Metadata,
    matrix: Matrix,
}
impl Profile {
    pub fn metadata(&self) -> &Metadata {
        &self.metadata
    }
    pub fn prepare(&self) -> Result<Prepared<'_>, Error> {
        Ok(Prepared { profile: self })
    }
    pub fn weight(&self, end: &str, partner: &str) -> Option<u64> {
        self.matrix.weight(end, partner)
    }
}
fn build_profile(i: usize, summary: Summary) -> Result<Arc<Profile>, Error> {
    let matrix = parse_published_matrix(SOURCES[i].raw, SOURCES[i].k)?;
    if matrix.source_blake3 != summary.provenance.source_blake3 {
        return Err(bad("embedded source hash changed"));
    }
    let data_revision =
        assembly_recipe::hash(&(MODEL, VERSION, &summary.profile, &matrix.weights))?;
    let metadata = Metadata {
        profile: summary.profile,
        version: summary.version,
        model: summary.model,
        name: summary.name,
        geometry_family: summary.geometry_family,
        enzyme_formulation: summary.enzyme_formulation,
        overhang_length: summary.overhang_length,
        conditions: summary.conditions,
        provenance: summary.provenance,
        source_sha256: summary.source_sha256,
        cells: summary.cells,
        data_revision,
    };
    Ok(Arc::new(Profile { metadata, matrix }))
}
pub fn load(reference: &ProfileRef) -> Result<Arc<Profile>, Error> {
    let (i, summary) = resolve(reference)?;
    PROFILES[i]
        .get_or_init(|| build_profile(i, summary).map_err(|e| e.to_string()))
        .as_ref()
        .cloned()
        .map_err(|e| bad(e))
}

pub fn source_bytes(reference: &ProfileRef) -> Result<&'static [u8], Error> {
    let (i, _) = resolve(reference)?;
    Ok(SOURCES[i].raw)
}
fn end_guard(k: usize, ends: &[String]) -> Result<(), Error> {
    if ends.is_empty()
        || ends.len() > numeric::MAX_ENDS
        || ends.iter().any(|s| index(s, k).is_none())
    {
        return Err(bad("choose1..256 concrete uppercase profile-sized ends"));
    }
    Ok(())
}
/// Score identity uses exact source/converter/conditions; it does not parse XLSX.
pub fn revision(reference: &ProfileRef, ends: &[String]) -> Result<String, Error> {
    let (_, s) = resolve(reference)?;
    end_guard(s.overhang_length, ends)?;
    revision_with_competitors(reference, ends, &[])
}
pub fn score_revision(reference: &ProfileRef, ends: &[String]) -> Result<String, Error> {
    revision(reference, ends)
}
pub fn revision_with_competitors(
    reference: &ProfileRef,
    joins: &[String],
    extras: &[String],
) -> Result<String, Error> {
    let (_, summary) = resolve(reference)?;
    numeric::physical_scope(summary.overhang_length, joins, extras)?;
    assembly_recipe::hash(&(
        MODEL,
        "literal-physical-competition",
        reference,
        joins,
        extras,
    ))
}
pub fn score_revision_with_competitors(
    reference: &ProfileRef,
    joins: &[String],
    extras: &[String],
) -> Result<String, Error> {
    revision_with_competitors(reference, joins, extras)
}
pub struct Prepared<'a> {
    profile: &'a Profile,
}
fn pairs(raw: Vec<numeric::PairScore>) -> Vec<PairScore> {
    raw.into_iter()
        .map(|p| PairScore {
            canonical_end: p.canonical_end,
            reverse_complement: p.reverse_complement,
            occurrences: p.occurrences,
            orientations: p
                .orientations
                .into_iter()
                .map(|r| Orientation {
                    end_5to3: r.end_5to3,
                    correct_partner_5to3: r.correct_partner_5to3,
                    correct_weight: r.correct_count,
                    total_weight: r.total_count,
                    missing_partners: r.missing_partners,
                })
                .collect(),
            estimate: match p.estimate {
                numeric::Estimate::Known {
                    correct_count,
                    total_count,
                    fidelity,
                } => Estimate::Known {
                    correct_weight: correct_count,
                    total_weight: total_count,
                    fidelity,
                },
                numeric::Estimate::Unknown { reason } => Estimate::Unknown { reason },
            },
        })
        .collect()
}
impl Prepared<'_> {
    pub fn score(&self, ends: &[String]) -> Result<Report, Error> {
        end_guard(self.profile.matrix.overhang_length, ends)?;
        self.score_with_competitors(ends, &[])
    }
    pub fn score_with_competitors(
        &self,
        ends: &[String],
        extras: &[String],
    ) -> Result<Report, Error> {
        let numeric = numeric::calculate_with_competitors(
            self.profile.matrix.overhang_length,
            ends,
            extras,
            |a, b| self.profile.weight(a, b),
        )?;
        Ok(Report{metadata:self.profile.metadata.clone(),context_revision:revision_with_competitors(&self.profile.metadata.profile,ends,extras)?,requested_ends:ends.to_vec(),physical_competitors:extras.to_vec(),physical_ends:numeric.physical_ends,pairs:pairs(numeric.pairs),overall_fidelity:numeric.overall_fidelity,log_overall_fidelity:numeric.log_overall_fidelity,
            warnings:vec!["Published symmetric ligation-frequency weights under the exact declared assay conditions; not raw independent directional sequencing counts, reaction yield, physical compatibility or experimental success.".into(),"Conditional pooled-physical-pair ranking. Each unique reverse-complement pair contributes once; both ordered rows and all present partners including self-mispairs enter its denominator. Extra actual5′→3′ ends are literal singleton competitors; their complements are not invented, and they are not multiplied as intended joins. Empty intended scope has no overall estimate. Equal partner weighting and independent pair outcomes are modeling assumptions; duplicates do not supply stoichiometry. Palindromes remain unknown.".into(),"All source table cells, including explicit zeros, are retained. No transposition, imputation, row normalization or pseudocount is introduced. A numerical1 is limited to these observed published weights; no uncertainty intervals or independent replicate model is inferred. Floating products can underflow; logarithms remain available for positive known ratios.".into(),"Data: CC-BY-4.0, Pryor et al.2020/PLOS; independent software: MIT. Exact enzyme formulation is not inferred from a shared recognition/cut geometry family.".into()]})
    }
}
pub fn score(reference: &ProfileRef, ends: &[String]) -> Result<Report, Error> {
    load(reference)?.prepare()?.score(ends)
}
pub fn score_with_competitors(
    reference: &ProfileRef,
    joins: &[String],
    extras: &[String],
) -> Result<Report, Error> {
    load(reference)?
        .prepare()?
        .score_with_competitors(joins, extras)
}
fn raw_pairs(scores: &[PairScore]) -> Result<Vec<numeric::PairScore>, Error> {
    if scores.len() > numeric::MAX_ENDS {
        return Err(Error::TooLarge);
    }
    scores
        .iter()
        .map(|p| {
            if p.canonical_end.len() > 4
                || p.reverse_complement.len() > 4
                || p.occurrences.len() > numeric::MAX_ENDS
                || p.orientations.len() > 2
            {
                return Err(Error::TooLarge);
            }
            let orientations = p
                .orientations
                .iter()
                .map(|r| {
                    if r.end_5to3.len() > 4
                        || r.correct_partner_5to3.len() > 4
                        || !r.missing_partners.is_empty()
                        || r.correct_weight.is_none()
                        || r.total_weight.is_none()
                    {
                        return Err(bad(
                            "complete published profile cannot have missing weight cells",
                        ));
                    }
                    Ok(numeric::Orientation {
                        end_5to3: r.end_5to3.clone(),
                        correct_partner_5to3: r.correct_partner_5to3.clone(),
                        correct_count: r.correct_weight,
                        total_count: r.total_weight,
                        missing_partners: vec![],
                    })
                })
                .collect::<Result<Vec<_>, Error>>()?;
            Ok(numeric::PairScore {
                canonical_end: p.canonical_end.clone(),
                reverse_complement: p.reverse_complement.clone(),
                occurrences: p.occurrences.clone(),
                orientations,
                estimate: match p.estimate {
                    Estimate::Known {
                        correct_weight,
                        total_weight,
                        fidelity,
                    } => numeric::Estimate::Known {
                        correct_count: correct_weight,
                        total_count: total_weight,
                        fidelity,
                    },
                    Estimate::Unknown { reason } => numeric::Estimate::Unknown { reason },
                },
            })
        })
        .collect()
}
impl Metadata {
    fn matches_summary(&self, s: &Summary) -> bool {
        self.profile == s.profile
            && self.version == s.version
            && self.model == s.model
            && self.name == s.name
            && self.geometry_family == s.geometry_family
            && self.enzyme_formulation == s.enzyme_formulation
            && self.overhang_length == s.overhang_length
            && self.conditions == s.conditions
            && self.provenance == s.provenance
            && self.source_sha256 == s.source_sha256
            && self.cells == s.cells
            && hex(&self.data_revision)
    }
}
impl Report {
    /// Cheap readonly source/converter/provenance and arithmetic consistency.
    /// Never parses XLSX, recreates a profile index or authenticates imported
    /// reported weights. Native scoring/ranking prepares the pinned source afresh.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if self.context_revision != expected
            || revision_with_competitors(
                &self.metadata.profile,
                &self.requested_ends,
                &self.physical_competitors,
            )? != expected
            || self.warnings.len() > 16
            || self.warnings.iter().any(|s| s.len() > 4096)
        {
            return Err(bad("score context/shape changed"));
        }
        let (_, summary) = resolve(&self.metadata.profile)?;
        if !self.metadata.matches_summary(&summary) {
            return Err(bad("source/version/conditions/license metadata changed"));
        }
        let raw = raw_pairs(&self.pairs)?;
        numeric::validate_numeric_with_competitors(
            self.metadata.overhang_length,
            &self.requested_ends,
            &self.physical_competitors,
            &self.physical_ends,
            &raw,
            self.overall_fidelity,
            self.log_overall_fidelity,
        )
    }
}
