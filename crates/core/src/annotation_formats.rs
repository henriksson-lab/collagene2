//! Independent, reference-bound GFF3/BED/GTF annotation exchange.
//! Original input is evidence; external annotations are never scientific verification.
use crate::{Alphabet, Error, Feature, Location, Molecule, Topology, formats::ExportReport};
use base64::{Engine, engine::general_purpose::STANDARD};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
mod gtf;

pub const MAX_ANNOTATION_BYTES: usize = 8 * 1024 * 1024;
pub const MAX_ANNOTATION_ROWS: usize = 10_000;
const MAX_PARTS: usize = 50_000;
const MAX_DIRECTIVES: usize = 20_000;
const MAX_REPORT_BYTES: usize = 32 * 1024 * 1024;
const EVIDENCE: &str = "annotation_exchange:source_record";
fn bad(message: impl Into<String>) -> Error {
    Error::Format(message.into())
}
fn unsupported(message: impl Into<String>) -> Error {
    Error::Unsupported(message.into())
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum AnnotationFormat {
    Gff3,
    Bed,
    Gtf,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct AnnotationProposal {
    pub id: String,
    pub feature: Feature,
    /// One-based physical input line numbers, including comments/directives.
    pub source_rows: Vec<usize>,
    pub warnings: Vec<String>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct AnnotationImportReport {
    pub schema_version: u16,
    pub reference_revision: String,
    pub format: AnnotationFormat,
    pub seqid: String,
    pub source_blake3: String,
    pub original_base64: String,
    pub proposals: Vec<AnnotationProposal>,
    pub skipped_seqids: Vec<String>,
    pub directives: Vec<String>,
    pub warnings: Vec<String>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Attribute {
    tag: String,
    values: Vec<String>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct GffRow {
    line: usize,
    seqid: String,
    source: String,
    kind: String,
    start: usize,
    end: usize,
    score: String,
    strand: String,
    phase: Option<u8>,
    attributes: Vec<Attribute>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct BedRow {
    line: usize,
    fields: Vec<String>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct SourceRecord {
    format: AnnotationFormat,
    kind: String,
    location: Location,
    gff: Vec<GffRow>,
    bed: Option<BedRow>,
    qualifier_hash: String,
}
fn number(s: &str) -> Result<usize, Error> {
    if s.is_empty() || !s.bytes().all(|b| b.is_ascii_digit()) {
        return Err(bad(
            "annotation coordinate is not an unsigned decimal integer",
        ));
    }
    s.parse()
        .map_err(|_| bad("annotation coordinate overflows this platform"))
}
fn text_token(s: &str, what: &str) -> Result<(), Error> {
    if s.is_empty() || s.len() > 8192 || s.chars().any(char::is_control) {
        return Err(bad(format!("invalid annotation {what}")));
    }
    Ok(())
}
fn seqid_token(s: &str) -> Result<(), Error> {
    text_token(s, "sequence ID")?;
    if s == "." {
        return Err(bad("an explicit non-missing sequence ID is required"));
    }
    Ok(())
}
fn hex(b: u8) -> Option<u8> {
    match b {
        b'0'..=b'9' => Some(b - b'0'),
        b'A'..=b'F' => Some(b - b'A' + 10),
        b'a'..=b'f' => Some(b - b'a' + 10),
        _ => None,
    }
}
fn decode(s: &str) -> Result<String, Error> {
    let mut out = Vec::with_capacity(s.len());
    let b = s.as_bytes();
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'%' {
            if i + 2 >= b.len() {
                return Err(bad("truncated GFF3 percent escape"));
            }
            out.push(
                hex(b[i + 1])
                    .zip(hex(b[i + 2]))
                    .map(|(a, b)| a * 16 + b)
                    .ok_or_else(|| bad("invalid GFF3 percent escape"))?,
            );
            i += 3;
        } else {
            out.push(b[i]);
            i += 1;
        }
    }
    String::from_utf8(out).map_err(|_| bad("GFF3 percent escapes are not valid UTF-8"))
}
fn encode(s: &str, seqid: bool) -> String {
    let mut out = String::new();
    for b in s.bytes() {
        let safe = if seqid {
            b.is_ascii_alphanumeric() || b".:^*$@!+_?-|".contains(&b)
        } else {
            b >= 32 && b < 127 && !b"%;=,&\t\r\n".contains(&b)
        };
        if safe {
            out.push(b as char)
        } else {
            out.push_str(&format!("%{b:02X}"));
        }
    }
    out
}
fn attributes(s: &str) -> Result<Vec<Attribute>, Error> {
    if s == "." {
        return Ok(vec![]);
    }
    let mut result = vec![];
    let mut tags = BTreeSet::new();
    for pair in s.split(';').filter(|s| !s.is_empty()) {
        let (tag, values) = pair.split_once('=').ok_or_else(|| {
            bad("GFF3 attributes require tag=value; GTF quotation syntax is unsupported")
        })?;
        if tag.contains([',', '&']) || values.contains(['=', '&']) {
            return Err(bad(
                "GFF3 reserved attribute characters must be percent-encoded",
            ));
        }
        let tag = decode(tag)?;
        text_token(&tag, "attribute tag")?;
        if !tags.insert(tag.clone()) {
            return Err(bad("duplicate GFF3 attribute tag"));
        }
        let values = values
            .split(',')
            .map(decode)
            .collect::<Result<Vec<_>, _>>()?;
        if values.len() > 1000 {
            return Err(Error::TooLarge);
        }
        if matches!(tag.as_str(), "ID" | "Name" | "Is_circular") && values.len() != 1 {
            return Err(bad("scalar GFF3 attribute has multiple values"));
        }
        if matches!(tag.as_str(), "ID" | "Parent" | "Derives_from") {
            for value in &values {
                text_token(value, "relationship identifier")?;
            }
        }
        result.push(Attribute { tag, values });
    }
    Ok(result)
}
fn attr<'a>(row: &'a GffRow, key: &str) -> Option<&'a [String]> {
    row.attributes
        .iter()
        .find(|a| a.tag == key)
        .map(|a| a.values.as_slice())
}
fn scalar(row: &GffRow, key: &str) -> Option<String> {
    attr(row, key).and_then(|v| v.first()).cloned()
}
fn location_parts(start: usize, end: usize, reference: &Molecule) -> Result<Vec<Location>, Error> {
    let n = reference.sequence.len();
    if start == 0 || end < start {
        return Err(bad("GFF3 coordinates require 1 <= start <= end"));
    }
    if start > n || end - start >= n {
        return Err(bad(
            "GFF3 location starts beyond the reference or spans more than one reference turn",
        ));
    }
    if end <= n {
        return Ok(vec![Location::range(start - 1, end)]);
    }
    if reference.topology != Topology::Circular {
        return Err(unsupported(
            "GFF3 end beyond reference requires an explicitly circular reference",
        ));
    }
    Ok(vec![
        Location::range(start - 1, n),
        Location::range(0, end - n),
    ])
}
fn join(parts: Vec<Location>) -> Location {
    if parts.len() == 1 {
        parts.into_iter().next().unwrap()
    } else {
        Location::Join(parts)
    }
}
fn check_input(reference: &Molecule, text: &str, seqid: &str) -> Result<(), Error> {
    reference.validate()?;
    seqid_token(seqid)?;
    if reference.sequence.alphabet() == Alphabet::Protein {
        return Err(unsupported(
            "annotation exchange requires a DNA/RNA reference",
        ));
    }
    if text.is_empty() {
        return Err(bad("annotation source is empty"));
    }
    if text.len() > MAX_ANNOTATION_BYTES {
        return Err(Error::TooLarge);
    }
    if text
        .bytes()
        .any(|b| b.is_ascii_control() && !matches!(b, b'\t' | b'\r' | b'\n'))
    {
        return Err(bad(
            "annotation source contains an unescaped control character",
        ));
    }
    if text.lines().any(|l| l.len() > 64 * 1024) {
        return Err(Error::TooLarge);
    }
    Ok(())
}
fn graph(rows: &[GffRow]) -> Result<(), Error> {
    let mut ids = BTreeSet::new();
    for row in rows {
        if let Some(id) = scalar(row, "ID") {
            ids.insert(id);
        }
    }
    let mut edges = BTreeMap::<String, BTreeSet<String>>::new();
    for row in rows {
        if let Some(parents) = attr(row, "Parent") {
            for parent in parents {
                if !ids.contains(parent) {
                    return Err(bad(format!("unresolved GFF3 Parent {parent:?}")));
                }
                if let Some(id) = scalar(row, "ID") {
                    edges.entry(id).or_default().insert(parent.clone());
                }
            }
        }
    }
    // Kahn traversal avoids recursive depth/stack growth and quadratic path enumeration.
    let mut remaining = ids
        .iter()
        .map(|id| (id.clone(), 0usize))
        .collect::<BTreeMap<_, _>>();
    let mut children = BTreeMap::<String, Vec<String>>::new();
    for (child, parents) in edges {
        *remaining.get_mut(&child).unwrap() = parents.len();
        for parent in parents {
            children.entry(parent).or_default().push(child.clone());
        }
    }
    let mut ready = remaining
        .iter()
        .filter(|(_, n)| **n == 0)
        .map(|(id, _)| id.clone())
        .collect::<Vec<_>>();
    let mut visited = 0;
    while let Some(id) = ready.pop() {
        visited += 1;
        for child in children.get(&id).into_iter().flatten() {
            let n = remaining.get_mut(child).unwrap();
            *n -= 1;
            if *n == 0 {
                ready.push(child.clone());
            }
        }
    }
    if visited != ids.len() {
        return Err(bad("cyclic GFF3 Parent relationships"));
    }
    Ok(())
}
fn parse_gff(
    reference: &Molecule,
    text: &str,
    seqid: &str,
) -> Result<(Vec<AnnotationProposal>, Vec<String>, Vec<String>), Error> {
    let mut rows = vec![];
    let mut directives = vec![];
    let mut seen_version = false;
    let mut fasta = None;
    let mut regions = BTreeMap::<String, (usize, usize)>::new();
    let mut seen_ids = BTreeSet::new();
    let mut completed = BTreeSet::new();
    let mut current_ids = BTreeSet::new();
    let mut since_boundary = 0usize;
    for (i, line) in text.lines().enumerate() {
        let line = line.strip_suffix('\r').unwrap_or(line);
        if line.is_empty() {
            continue;
        }
        if line == "##FASTA" {
            directives.push(line.into());
            fasta = Some(text.lines().skip(i + 1).collect::<Vec<_>>().join("\n"));
            break;
        }
        if line.starts_with('#') {
            if directives.len() == MAX_DIRECTIVES {
                return Err(Error::TooLarge);
            }
            directives.push(line.into());
            if line.starts_with("##gff-version") {
                if line != "##gff-version 3" || seen_version || !rows.is_empty() {
                    return Err(unsupported(
                        "GFF3 requires one leading ##gff-version 3 directive",
                    ));
                }
                seen_version = true;
            }
            if let Some(region) = line.strip_prefix("##sequence-region ") {
                let fields = region.split_ascii_whitespace().collect::<Vec<_>>();
                if fields.len() != 3 {
                    return Err(bad("invalid GFF3 sequence-region directive"));
                }
                let id = decode(fields[0])?;
                seqid_token(&id)?;
                let (a, b) = (number(fields[1])?, number(fields[2])?);
                if a == 0 || a > b || regions.insert(id.clone(), (a, b)).is_some() {
                    return Err(bad("invalid or duplicate GFF3 sequence-region"));
                }
                if id == seqid && b > reference.sequence.len() {
                    return Err(bad("GFF3 sequence-region exceeds mapped reference"));
                }
            }
            if line == "###" {
                for row in &rows[since_boundary..] {
                    let row: &GffRow = row;
                    for parent in attr(row, "Parent").into_iter().flatten() {
                        if !seen_ids.contains(parent) {
                            return Err(bad(
                                "GFF3 ### occurs before a referenced Parent is defined",
                            ));
                        }
                    }
                }
                completed.extend(current_ids.iter().cloned());
                current_ids.clear();
                since_boundary = rows.len();
            }
            continue;
        }
        if !seen_version {
            return Err(bad("missing leading GFF3 version directive"));
        }
        if rows.len() == MAX_ANNOTATION_ROWS {
            return Err(Error::TooLarge);
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.len() != 9 {
            return Err(bad(format!(
                "GFF3 line {} must have nine tab-separated fields",
                i + 1
            )));
        }
        if fields[0]
            .bytes()
            .any(|b| !(b.is_ascii_alphanumeric() || b".:^*$@!+_?-|%".contains(&b)))
        {
            return Err(bad(
                "GFF3 sequence ID contains an unescaped reserved character",
            ));
        }
        let id = decode(fields[0])?;
        seqid_token(&id)?;
        let (start, end) = (number(fields[3])?, number(fields[4])?);
        if start == 0 || end < start {
            return Err(bad("GFF3 start/end are not positive ordered coordinates"));
        }
        if fields[5] != "." && fields[5].parse::<f64>().ok().is_none_or(|n| !n.is_finite()) {
            return Err(bad("GFF3 score must be a finite number or ."));
        }
        if !matches!(fields[6], "+" | "-" | "." | "?") {
            return Err(bad("invalid GFF3 strand"));
        }
        let phase = match fields[7] {
            "." => None,
            "0" => Some(0),
            "1" => Some(1),
            "2" => Some(2),
            _ => return Err(bad("GFF3 phase must be 0, 1, 2 or .")),
        };
        let kind = decode(fields[2])?;
        text_token(&kind, "type")?;
        if kind == "." {
            return Err(bad("GFF3 feature type is required"));
        }
        if kind == "CDS" && phase.is_none() {
            return Err(bad("GFF3 CDS requires an explicit phase"));
        }
        if kind != "CDS" && phase.is_some() {
            return Err(bad("GFF3 phase is defined only for CDS"));
        }
        let row = GffRow {
            line: i + 1,
            seqid: id,
            source: decode(fields[1])?,
            kind,
            start,
            end,
            score: fields[5].into(),
            strand: fields[6].into(),
            phase,
            attributes: attributes(fields[8])?,
        };
        if let Some(id) = scalar(&row, "ID") {
            if completed.contains(&id) {
                return Err(bad(
                    "GFF3 ID reappears after a ### completed-feature boundary",
                ));
            }
            current_ids.insert(id.clone());
            seen_ids.insert(id);
        }
        if let Some(&(a, b)) = regions.get(&row.seqid) {
            if start < a
                || (end > b
                    && !(row.seqid == seqid
                        && reference.topology == Topology::Circular
                        && end - b <= reference.sequence.len()))
            {
                return Err(bad("GFF3 row exceeds its declared sequence-region"));
            }
        }
        if row.seqid == seqid {
            location_parts(start, end, reference)?;
            if scalar(&row, "Is_circular").as_deref() == Some("true")
                && reference.topology != Topology::Circular
            {
                return Err(bad(
                    "GFF3 circular landmark conflicts with the mapped linear reference",
                ));
            }
        }
        rows.push(row);
    }
    if !seen_version {
        return Err(bad("missing GFF3 version directive"));
    }
    if let Some(fasta) = fasta {
        let report = crate::formats::import_text(&fasta)?;
        let mut names = BTreeSet::new();
        for molecule in report.molecules {
            let id = decode(&molecule.name)?;
            if !names.insert(id.clone()) {
                return Err(bad("duplicate embedded GFF3 FASTA identifier"));
            }
            if id == seqid && molecule.sequence != reference.sequence {
                return Err(bad(
                    "embedded GFF3 FASTA disagrees with mapped reference sequence/alphabet",
                ));
            }
        }
    }
    graph(&rows)?;
    let skipped = rows
        .iter()
        .filter(|r| r.seqid != seqid)
        .map(|r| r.seqid.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let mut groups = Vec::<Vec<GffRow>>::new();
    let mut group_ids = BTreeMap::<String, usize>::new();
    for row in rows {
        if let Some(id) = scalar(&row, "ID") {
            if let Some(&index) = group_ids.get(&id) {
                groups[index].push(row);
            } else {
                group_ids.insert(id, groups.len());
                groups.push(vec![row]);
            }
        } else {
            groups.push(vec![row]);
        }
    }
    let mut proposals = vec![];
    for group in groups {
        if !group.iter().any(|r| r.seqid == seqid) {
            continue;
        }
        proposals.push(gff_proposal(group, reference)?);
    }
    Ok((proposals, skipped, directives))
}

fn gff_proposal(mut group: Vec<GffRow>, reference: &Molecule) -> Result<AnnotationProposal, Error> {
    if group.is_empty() || group.len() > MAX_ANNOTATION_ROWS {
        return Err(bad("invalid retained GFF3 row group"));
    }
    for row in &group {
        if !matches!(row.strand.as_str(), "+" | "-" | "." | "?")
            || row.phase.is_some_and(|p| p > 2)
            || (row.kind == "CDS" && row.phase.is_none())
            || (row.kind != "CDS" && row.phase.is_some())
        {
            return Err(bad("invalid retained GFF3 strand/phase"));
        }
    }
    if group.len() > 1 {
        let id =
            scalar(&group[0], "ID").ok_or_else(|| bad("discontinuous GFF3 feature lacks ID"))?;
        if group
            .iter()
            .any(|r| scalar(r, "ID").as_deref() != Some(id.as_str()))
        {
            return Err(bad("discontinuous GFF3 group has different IDs"));
        }
    }
    group.sort_by_key(|r| (r.start, r.end, r.line));
    let first = group[0].clone();
    if group
        .iter()
        .any(|r| r.seqid != first.seqid || r.kind != first.kind || r.strand != first.strand)
    {
        return Err(unsupported(
            "one GFF3 ID spans different references, types or strands; explicit modeling is required",
        ));
    }
    let mut warnings = vec![
        "Imported annotation requires source/biological review; this is not sequence verification."
            .into(),
    ];
    let mut parts = vec![];
    let mut previous = 0;
    for row in &group {
        if row.start <= previous {
            return Err(bad(
                "overlapping/disordered parts of one GFF3 ID are not a unique compound annotation",
            ));
        }
        previous = row.end;
        parts.extend(location_parts(row.start, row.end, reference)?);
    }
    let interval_kind = matches!(
        first.kind.as_str(),
        "gene"
            | "pseudogene"
            | "CDS"
            | "exon"
            | "intron"
            | "mRNA"
            | "transcript"
            | "promoter"
            | "terminator"
            | "sequence_feature"
            | "SNP"
            | "SNV"
            | "single_nucleotide_variant"
    );
    let point = !interval_kind && group.iter().any(|r| r.start == r.end);
    let wrapped = group.iter().any(|r| r.end > reference.sequence.len());
    let circular_anchor = reference.topology == Topology::Circular
        && group.len() > 1
        && (wrapped
            || (group[0].start == 1 && group.last().unwrap().end == reference.sequence.len()));
    let location = if point {
        warnings.push("GFF3 equal start/end with a point-site or unrecognized type has unresolved interval/boundary semantics. Known interval types use ordinary one-based inclusive single-base coordinates. Original rows remain available for review.".into());
        Location::Unresolved(format!(
            "GFF3 point feature at source lines {:?}",
            group.iter().map(|r| r.line).collect::<Vec<_>>()
        ))
    } else if circular_anchor {
        warnings.push("Multi-row circular origin anchor is unresolved; original parts/phase remain source evidence.".into());
        Location::Unresolved("GFF3 circular compound biological anchor requires review".into())
    } else if matches!(first.strand.as_str(), "." | "?") {
        warnings.push("Strand is unstranded/unknown. Reference coordinates are retained, but exact biological extraction is disabled.".into());
        Location::Order(parts)
    } else {
        let location = join(parts);
        if first.strand == "-" {
            Location::Complement(Box::new(location))
        } else {
            location
        }
    };
    let mut qualifiers = vec![
        (
            "label".into(),
            Some(
                scalar(&first, "Name")
                    .or_else(|| scalar(&first, "ID"))
                    .unwrap_or_else(|| first.kind.clone()),
            ),
        ),
        ("gff3:strand".into(), Some(first.strand.clone())),
    ];
    if let Some(id) = scalar(&first, "ID") {
        qualifiers.push(("gff3:ID".into(), Some(id)));
    }
    for row in &group {
        for a in &row.attributes {
            for v in &a.values {
                let q = (format!("gff3:attribute:{}", a.tag), Some(v.clone()));
                if !qualifiers.contains(&q) {
                    qualifiers.push(q);
                }
            }
        }
    }
    if group
        .iter()
        .any(|r| attr(r, "Target").is_some() || attr(r, "Gap").is_some())
    {
        warnings.push(
            "Target/Gap are retained alignment evidence, not a reconstructed pairwise alignment."
                .into(),
        );
    }
    if first.kind == "CDS"
        && !point
        && !circular_anchor
        && matches!(first.strand.as_str(), "+" | "-")
    {
        let biological = if first.strand == "-" {
            group.iter().rev().collect::<Vec<_>>()
        } else {
            group.iter().collect()
        };
        let initial = biological[0].phase.unwrap();
        let mut length = 0usize;
        let mut consistent = true;
        for (i, row) in biological.iter().enumerate() {
            if i > 0 {
                let remainder = (length % 3 + 3 - initial as usize) % 3;
                let expected = ((3 - remainder) % 3) as u8;
                if row.phase != Some(expected) {
                    consistent = false;
                }
            }
            length += row.end - row.start + 1;
        }
        if consistent {
            qualifiers.push(("codon_start".into(), Some((initial + 1).to_string())));
        } else {
            warnings.push("Discontinuous CDS phases are inconsistent; no current codon_start or translation is asserted.".into());
            qualifiers.push(("gff3:inconsistent_phase".into(), Some("true".into())));
        }
    }
    let evidence = SourceRecord {
        format: AnnotationFormat::Gff3,
        kind: first.kind.clone(),
        location: location.clone(),
        gff: group.clone(),
        bed: None,
        qualifier_hash: qualifier_hash(&qualifiers)?,
    };
    qualifiers.push((
        EVIDENCE.into(),
        Some(serde_json::to_string(&evidence).map_err(|e| bad(e.to_string()))?),
    ));
    Ok(AnnotationProposal {
        id: String::new(),
        feature: Feature {
            kind: first.kind,
            location,
            qualifiers,
            needs_review: true,
        },
        source_rows: group.iter().map(|r| r.line).collect(),
        warnings,
    })
}

fn list(s: &str) -> Result<Vec<usize>, Error> {
    let s = s.strip_suffix(',').unwrap_or(s);
    if s.is_empty() {
        return Err(bad("empty BED block list"));
    }
    s.split(',').map(number).collect()
}
fn parse_bed(
    reference: &Molecule,
    text: &str,
    seqid: &str,
) -> Result<(Vec<AnnotationProposal>, Vec<String>, Vec<String>), Error> {
    let mut proposals = vec![];
    let mut skipped = BTreeSet::new();
    let mut directives = vec![];
    let mut width = None;
    let mut count = 0;
    let mut parts_count = 0;
    for (i, line) in text.lines().enumerate() {
        let line = line.strip_suffix('\r').unwrap_or(line);
        if line.is_empty() {
            continue;
        }
        if line.starts_with('#')
            || line == "track"
            || line.starts_with("track ")
            || line.starts_with("track\t")
            || line == "browser"
            || line.starts_with("browser ")
            || line.starts_with("browser\t")
        {
            if line.starts_with("track") {
                width = None;
            }
            if directives.len() == MAX_DIRECTIVES {
                return Err(Error::TooLarge);
            }
            directives.push(line.into());
            continue;
        }
        count += 1;
        if count > MAX_ANNOTATION_ROWS {
            return Err(Error::TooLarge);
        }
        let f = if line.contains('\t') {
            line.split('\t').collect::<Vec<_>>()
        } else {
            line.split_ascii_whitespace().collect::<Vec<_>>()
        };
        if !(3..=9).contains(&f.len()) && f.len() != 12 {
            return Err(unsupported(
                "BED requires 3–9 or all 12 fields; partial block columns and BED extensions are unsupported",
            ));
        }
        if width.is_some_and(|n| n != f.len()) {
            return Err(bad("BED field count changes within a track"));
        }
        width = Some(f.len());
        seqid_token(f[0])?;
        let (a, b) = (number(f[1])?, number(f[2])?);
        if a > b {
            return Err(bad("BED requires chromStart <= chromEnd"));
        }
        if f.len() >= 5 && number(f[4])? > 1000 {
            return Err(bad("BED score exceeds 1000"));
        }
        let strand = if f.len() >= 6 { f[5] } else { "." };
        if !matches!(strand, "+" | "-" | ".") {
            return Err(bad("invalid BED strand"));
        }
        if f.len() >= 7 {
            let start = number(f[6])?;
            if start < a || start > b {
                return Err(bad("BED thickStart outside interval"));
            }
            if f.len() >= 8 {
                let end = number(f[7])?;
                if end < start || end > b {
                    return Err(bad("BED thickEnd outside interval"));
                }
            }
        }
        if f.len() >= 9 && f[8] != "0" {
            let rgb = list(f[8])?;
            if rgb.len() != 3 || rgb.iter().any(|n| *n > 255) {
                return Err(bad("invalid BED RGB color"));
            }
        }
        let mut spans = vec![(a, b)];
        if f.len() == 12 {
            let n = number(f[9])?;
            if n == 0 || n > MAX_PARTS {
                return Err(bad("invalid BED block count"));
            }
            let sizes = list(f[10])?;
            let starts = list(f[11])?;
            if sizes.len() != n || starts.len() != n || starts[0] != 0 {
                return Err(bad(
                    "BED block sizes/starts/count disagree or first block does not start at chromStart",
                ));
            }
            spans.clear();
            let mut previous = 0;
            for (j, (&s, &l)) in starts.iter().zip(&sizes).enumerate() {
                let end = s.checked_add(l).ok_or(Error::TooLarge)?;
                if l == 0 || (j > 0 && s < previous) || end > b - a {
                    return Err(bad("BED blocks overlap, are empty or exceed interval"));
                }
                spans.push((
                    a.checked_add(s).ok_or(Error::TooLarge)?,
                    a.checked_add(end).ok_or(Error::TooLarge)?,
                ));
                previous = end;
            }
            if previous != b - a {
                return Err(bad("last BED block does not end at chromEnd"));
            }
        }
        parts_count += spans.len();
        if parts_count > MAX_PARTS {
            return Err(Error::TooLarge);
        }
        if f[0] != seqid {
            skipped.insert(f[0].to_string());
            continue;
        }
        if b > reference.sequence.len() {
            return Err(bad("BED interval exceeds explicitly mapped reference"));
        }
        let mut warnings=vec!["Imported annotation requires source/biological review; BED thick regions do not establish CDS or phase.".into()];
        let parts = spans
            .into_iter()
            .map(|(a, b)| Location::range(a, b))
            .collect::<Vec<_>>();
        let location = if a == b {
            warnings.push("BED zero-length insertion boundary retained unresolved; no base interval is invented.".into());
            Location::Unresolved(format!("BED boundary {a} on {seqid}"))
        } else if strand == "." {
            warnings.push("BED strand is not specified; coordinates are retained while biological extraction is disabled.".into());
            Location::Order(parts)
        } else {
            let location = join(parts);
            if strand == "-" {
                Location::Complement(Box::new(location))
            } else {
                location
            }
        };
        let mut qualifiers = vec![
            (
                "label".into(),
                Some(if f.len() >= 4 {
                    f[3].into()
                } else {
                    format!("{seqid}:{a}-{b}")
                }),
            ),
            ("bed:strand".into(), Some(strand.into())),
        ];
        let evidence = SourceRecord {
            format: AnnotationFormat::Bed,
            kind: "misc_feature".into(),
            location: location.clone(),
            gff: vec![],
            bed: Some(BedRow {
                line: i + 1,
                fields: f.iter().map(|s| s.to_string()).collect(),
            }),
            qualifier_hash: qualifier_hash(&qualifiers)?,
        };
        qualifiers.push((
            EVIDENCE.into(),
            Some(serde_json::to_string(&evidence).map_err(|e| bad(e.to_string()))?),
        ));
        proposals.push(AnnotationProposal {
            id: String::new(),
            feature: Feature {
                kind: "misc_feature".into(),
                location,
                qualifiers,
                needs_review: true,
            },
            source_rows: vec![i + 1],
            warnings,
        });
    }
    Ok((proposals, skipped.into_iter().collect(), directives))
}
struct Limit(usize);
impl std::io::Write for Limit {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.0 = self
            .0
            .checked_add(b.len())
            .filter(|n| *n <= MAX_REPORT_BYTES)
            .ok_or_else(|| std::io::Error::other("annotation output exceeds limit"))?;
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn bounded(value: &impl Serialize) -> Result<(), Error> {
    serde_json::to_writer(&mut Limit(0), value).map_err(|_| Error::TooLarge)
}
pub fn import_annotations(
    reference: &Molecule,
    text: &str,
    format: AnnotationFormat,
    seqid: &str,
) -> Result<AnnotationImportReport, Error> {
    check_input(reference, text, seqid)?;
    let (mut proposals, skipped_seqids, directives) = match format {
        AnnotationFormat::Gff3 => parse_gff(reference, text, seqid)?,
        AnnotationFormat::Bed => parse_bed(reference, text, seqid)?,
        AnnotationFormat::Gtf => gtf::parse(reference, text, seqid)?,
    };
    if proposals.is_empty() {
        return Err(bad(
            "annotation file has no rows for the explicitly selected sequence ID",
        ));
    }
    let hash = blake3::hash(text.as_bytes()).to_hex().to_string();
    for p in &mut proposals {
        p.feature.location.validate(reference.sequence.len())?;
        let tuple = (&hash, format, seqid, &p.source_rows, &p.feature);
        p.id = format!(
            "annotation-{}",
            blake3::hash(&serde_json::to_vec(&tuple).map_err(|e| bad(e.to_string()))?).to_hex()
        );
    }
    let mut warnings=vec!["Source sequence ID is an explicit coordinate mapping chosen by the caller; annotations do not establish reference identity or experimental verification.".into(),"Raw annotation source and row evidence are retained; unknown ontology terms/relationships and visualization directives are not executed.".into()];
    if !skipped_seqids.is_empty() {
        warnings.push(format!("{} other sequence IDs remain in original-source evidence and are not applied to this reference.",skipped_seqids.len()));
    }
    let report = AnnotationImportReport {
        schema_version: 1,
        reference_revision: reference.revision(),
        format,
        seqid: seqid.into(),
        source_blake3: hash,
        original_base64: STANDARD.encode(text),
        proposals,
        skipped_seqids,
        directives,
        warnings,
    };
    bounded(&report)?;
    Ok(report)
}
impl AnnotationImportReport {
    pub fn original_bytes(&self) -> Result<Vec<u8>, Error> {
        if self.original_base64.len() > MAX_ANNOTATION_BYTES.div_ceil(3) * 4 {
            return Err(Error::TooLarge);
        }
        let bytes = STANDARD
            .decode(&self.original_base64)
            .map_err(|_| bad("invalid annotation original-source base64"))?;
        if bytes.len() > MAX_ANNOTATION_BYTES
            || blake3::hash(&bytes).to_hex().as_str() != self.source_blake3
        {
            return Err(bad("annotation original-source hash mismatch"));
        }
        Ok(bytes)
    }
    pub fn validate(&self, reference: &Molecule) -> Result<(), Error> {
        bounded(self)?;
        if self.schema_version != 1 {
            return Err(unsupported("unknown annotation import schema"));
        }
        if self.reference_revision != reference.revision() {
            return Err(bad(
                "annotation import belongs to a different reference revision",
            ));
        }
        let raw = self.original_bytes()?;
        let text = std::str::from_utf8(&raw).map_err(|_| bad("annotation source is not UTF-8"))?;
        let rebuilt = import_annotations(reference, text, self.format, &self.seqid)?;
        if self != &rebuilt {
            return Err(bad(
                "annotation proposals/evidence disagree with retained source and reference mapping",
            ));
        }
        Ok(())
    }
}
pub fn apply_annotations(
    reference: &Molecule,
    report: &AnnotationImportReport,
    selection: &[usize],
) -> Result<Molecule, Error> {
    report.validate(reference)?;
    if selection.is_empty() {
        return Err(bad("select at least one annotation proposal"));
    }
    if selection.len() > MAX_ANNOTATION_ROWS {
        return Err(Error::TooLarge);
    }
    let mut seen = BTreeSet::new();
    let mut additions = vec![];
    for &i in selection {
        if !seen.insert(i) {
            return Err(bad("duplicate annotation selection"));
        }
        let p = report
            .proposals
            .get(i)
            .ok_or_else(|| bad("annotation selection is out of range"))?;
        if reference.features.contains(&p.feature) {
            return Err(bad(
                "selected annotation is already present; no duplicate was added",
            ));
        }
        additions.push(p.feature.clone());
    }
    let mut result = reference.clone();
    result.features.extend(additions);
    let prefix = format!("annotation_import:{}", report.source_blake3);
    if let Some(existing) = result.metadata.get(&format!("{prefix}:original_base64")) {
        if existing != &report.original_base64 {
            return Err(bad("conflicting annotation source evidence"));
        }
    }
    result.metadata.insert(
        format!("{prefix}:original_base64"),
        report.original_base64.clone(),
    );
    result
        .metadata
        .insert(format!("{prefix}:format"), format!("{:?}", report.format));
    let key = format!("{prefix}:applications");
    let mut applications = result
        .metadata
        .get(&key)
        .map(|s| {
            serde_json::from_str::<Vec<(String, String, Vec<String>)>>(s)
                .map_err(|_| bad("invalid previous annotation import provenance"))
        })
        .transpose()?
        .unwrap_or_default();
    applications.push((
        report.reference_revision.clone(),
        report.seqid.clone(),
        selection
            .iter()
            .map(|&i| report.proposals[i].id.clone())
            .collect(),
    ));
    result.metadata.insert(
        key,
        serde_json::to_string(&applications).map_err(|e| bad(e.to_string()))?,
    );
    result.validate()?;
    crate::protocol::save_document(&result)?;
    Ok(result)
}

fn qualifier_hash(qualifiers: &[(String, Option<String>)]) -> Result<String, Error> {
    let values = qualifiers
        .iter()
        .filter(|(k, _)| k != EVIDENCE)
        .collect::<Vec<_>>();
    Ok(
        blake3::hash(&serde_json::to_vec(&values).map_err(|e| bad(e.to_string()))?)
            .to_hex()
            .to_string(),
    )
}
fn exact_parts(location: &Location) -> Result<(Vec<(usize, usize)>, bool), Error> {
    let mut inner = location;
    let mut reverse = false;
    while let Location::Complement(next) = inner {
        reverse = !reverse;
        inner = next;
    }
    let candidates = match inner {
        Location::Join(parts) => parts.iter().collect::<Vec<_>>(),
        Location::Range { .. } => vec![inner],
        _ => {
            return Err(unsupported(
                "annotation export requires exact ordered intervals; unresolved/unknown-strand/order locations need review",
            ));
        }
    };
    let mut result = vec![];
    for candidate in candidates {
        match candidate {
            Location::Range {
                start,
                end,
                before: false,
                after: false,
            } => result.push((*start, *end)),
            _ => {
                return Err(unsupported(
                    "mixed-strand or fuzzy annotation components cannot be exported as one GFF3/BED feature",
                ));
            }
        }
    }
    if result.len() > MAX_PARTS {
        return Err(Error::TooLarge);
    }
    Ok((result, reverse))
}
fn push(out: &mut String, text: &str) -> Result<(), Error> {
    if out
        .len()
        .checked_add(text.len())
        .is_none_or(|n| n > MAX_REPORT_BYTES)
    {
        return Err(Error::TooLarge);
    }
    out.push_str(text);
    Ok(())
}
fn gff_line(
    seqid: &str,
    source: &str,
    kind: &str,
    start: usize,
    end: usize,
    score: &str,
    strand: &str,
    phase: Option<u8>,
    attributes: &[Attribute],
) -> String {
    let a = if attributes.is_empty() {
        ".".into()
    } else {
        attributes
            .iter()
            .map(|a| {
                format!(
                    "{}={}",
                    encode(&a.tag, false),
                    a.values
                        .iter()
                        .map(|v| encode(v, false))
                        .collect::<Vec<_>>()
                        .join(",")
                )
            })
            .collect::<Vec<_>>()
            .join(";")
    };
    format!(
        "{}\t{}\t{}\t{start}\t{end}\t{score}\t{strand}\t{}\t{a}\n",
        encode(seqid, true),
        encode(source, false),
        encode(kind, false),
        phase.map(|p| p.to_string()).unwrap_or_else(|| ".".into())
    )
}
/// Export current supported features. Original-source downloads are a separate operation.
pub fn export_annotations(
    reference: &Molecule,
    format: AnnotationFormat,
    seqid: &str,
) -> Result<ExportReport, Error> {
    reference.validate()?;
    seqid_token(seqid)?;
    if format == AnnotationFormat::Bed && seqid.chars().any(char::is_whitespace) {
        return Err(bad("BED sequence ID cannot contain whitespace"));
    }

    if reference.sequence.alphabet() == Alphabet::Protein {
        return Err(unsupported(
            "annotation exchange requires a nucleotide reference",
        ));
    }
    if reference.features.len() > MAX_ANNOTATION_ROWS {
        return Err(Error::TooLarge);
    }
    let mut out = String::new();
    let mut warnings=vec!["Annotation exports do not include the reference sequence, original attachments or scientific verification status; preserve the native workspace for full evidence.".into()];
    warnings.push("Original track/browser directives remain original-source evidence and are not replayed as current visualization settings.".into());
    if format == AnnotationFormat::Gff3 {
        push(
            &mut out,
            &format!(
                "##gff-version 3\n##sequence-region {} 1 {}\n",
                encode(seqid, true),
                reference.sequence.len()
            ),
        )?;
    }
    // BED output is always BED12 to avoid mixing row widths. Source BED3–9 rows
    // are normalized using explicit visualization placeholders, never CDS inference.
    for (index, feature) in reference.features.iter().enumerate() {
        let evidence = feature
            .qualifiers
            .iter()
            .filter(|(k, _)| k == EVIDENCE)
            .collect::<Vec<_>>();
        if evidence.len() > 1 {
            return Err(bad("duplicate annotation source-record evidence"));
        }
        let source = if let Some((_, Some(raw))) = evidence.first().copied() {
            let source: SourceRecord = serde_json::from_str(raw)
                .map_err(|_| bad("invalid retained annotation source record"))?;
            if source.kind != feature.kind
                || source.location != feature.location
                || source.qualifier_hash != qualifier_hash(&feature.qualifiers)?
            {
                return Err(unsupported(
                    "imported annotation changed since its source record: current phase/relationship/visualization remapping requires review; original source remains available",
                ));
            }
            let reconstructed = match source.format {
                AnnotationFormat::Gff3 => gff_proposal(source.gff.clone(), reference)?.feature,
                AnnotationFormat::Gtf => {
                    if source.gff.len() != 1 || source.bed.is_some() {
                        return Err(bad("GTF evidence must retain exactly one row"));
                    }
                    gtf::proposal(source.gff[0].clone(), reference)?.feature
                }
                AnnotationFormat::Bed => {
                    let row = source
                        .bed
                        .as_ref()
                        .ok_or_else(|| bad("missing retained BED record"))?;
                    let id = row
                        .fields
                        .first()
                        .ok_or_else(|| bad("empty retained BED record"))?;
                    let parsed = parse_bed(reference, &row.fields.join("\t"), id)?;
                    if parsed.0.len() != 1 {
                        return Err(bad("invalid retained BED record"));
                    }
                    parsed.0.into_iter().next().unwrap().feature
                }
            };
            if reconstructed.kind != feature.kind
                || reconstructed.location != feature.location
                || qualifier_hash(&reconstructed.qualifiers)?
                    != qualifier_hash(&feature.qualifiers)?
            {
                return Err(bad(
                    "retained annotation rows disagree with current feature geometry/phase/qualifiers",
                ));
            }
            Some(source)
        } else {
            None
        };
        if let Some(source) = &source {
            if source.format == AnnotationFormat::Gtf && format == AnnotationFormat::Gtf {
                push(&mut out, &gtf::line(seqid, &source.gff[0])?)?;
                warnings.push(format!("Feature {} preserves GTF row coordinates, strand, frame and repeated attributes as source evidence. Quotation/numeric lexical form is normalized; no transcript or splice model is inferred.",index+1));
                continue;
            }
        }
        if format == AnnotationFormat::Gtf {
            gtf::generic(feature, reference, seqid, index, &mut out, &mut warnings)?;
            continue;
        }
        if let Some(source) = &source {
            if source.format == AnnotationFormat::Gff3 && format == AnnotationFormat::Gff3 {
                for row in &source.gff {
                    push(
                        &mut out,
                        &gff_line(
                            seqid,
                            &row.source,
                            &row.kind,
                            row.start,
                            row.end,
                            &row.score,
                            &row.strand,
                            row.phase,
                            &row.attributes,
                        ),
                    )?;
                }
                warnings.push(format!("Feature {} preserves recorded GFF3 row phases, attributes and relationships as source evidence; these are not independently revalidated biological assertions.",index+1));
                continue;
            }
        }
        if let Some(source) = &source {
            if source.format == AnnotationFormat::Bed && format == AnnotationFormat::Bed {
                let row = source
                    .bed
                    .as_ref()
                    .ok_or_else(|| bad("retained BED source has no row"))?;
                let mut f = row.fields.clone();
                f[0] = seqid.into();
                if f.len() == 12 {
                    push(&mut out, &format!("{}\n", f.join("\t")))?;
                    continue;
                }
                let a = number(&f[1])?;
                let b = number(&f[2])?;
                if a == b {
                    return Err(unsupported(
                        "zero-length BED source cannot be normalized into the BED12 output profile; download original source",
                    ));
                }
                let name = f.get(3).cloned().unwrap_or_else(|| feature.label().into());
                let score = f.get(4).cloned().unwrap_or_else(|| "0".into());
                let strand = f.get(5).cloned().unwrap_or_else(|| ".".into());
                let thick_start = f.get(6).cloned().unwrap_or_else(|| a.to_string());
                let thick_end = f.get(7).cloned().unwrap_or_else(|| thick_start.clone());
                let color = f.get(8).cloned().unwrap_or_else(|| "0".into());
                push(
                    &mut out,
                    &format!(
                        "{seqid}\t{a}\t{b}\t{name}\t{score}\t{strand}\t{thick_start}\t{thick_end}\t{color}\t1\t{}\t0\n",
                        b - a
                    ),
                )?;
                warnings.push(format!("Feature {} normalized BED{} to BED12. Missing score/color/thick fields use 0/no-thick display placeholders, not measured score or inferred coding sequence.",index+1,f.len()));
                continue;
            }
        }
        let (mut parts, reverse) = exact_parts(&feature.location)?;
        let origin = parts.len() == 2
            && parts[0].1 == reference.sequence.len()
            && parts[1].0 == 0
            && parts[0].0 > parts[1].1;
        if parts.windows(2).any(|w| w[0].1 > w[1].0)
            && !(format == AnnotationFormat::Gff3
                && origin
                && reference.topology == Topology::Circular)
        {
            return Err(unsupported(
                "overlapping/non-reference-ordered or origin-crossing annotation cannot be flattened into this format",
            ));
        }
        if format == AnnotationFormat::Bed {
            let a = parts[0].0;
            let b = parts.last().unwrap().1;
            let name = feature.label();
            text_token(name, "BED label")?;
            if name.contains(['\t', '\r', '\n']) {
                return Err(bad("BED labels cannot contain field separators"));
            }
            let sizes = parts
                .iter()
                .map(|(a, b)| (b - a).to_string())
                .collect::<Vec<_>>()
                .join(",");
            let starts = parts
                .iter()
                .map(|(start, _)| (start - a).to_string())
                .collect::<Vec<_>>()
                .join(",");
            push(
                &mut out,
                &format!(
                    "{seqid}\t{a}\t{b}\t{name}\t0\t{}\t{a}\t{a}\t0\t{}\t{sizes}\t{starts}\n",
                    if reverse { "-" } else { "+" },
                    parts.len()
                ),
            )?;
            warnings.push(format!("Feature {} exported as BED blocks: feature type, biological qualifiers and CDS phase are omitted; score/color 0 and no-thick extent are display placeholders.",index+1));
        } else {
            let kind = if matches!(
                feature.kind.as_str(),
                "gene"
                    | "pseudogene"
                    | "CDS"
                    | "mRNA"
                    | "transcript"
                    | "exon"
                    | "intron"
                    | "promoter"
                    | "terminator"
                    | "sequence_feature"
                    | "SNP"
                    | "SNV"
                    | "single_nucleotide_variant"
            ) {
                feature.kind.as_str()
            } else {
                "sequence_feature"
            };
            let id = feature
                .qualifiers
                .iter()
                .find(|(k, _)| k == "gff3:ID")
                .and_then(|(_, v)| v.clone())
                .unwrap_or_else(|| format!("cg2-feature-{}", index + 1));
            let mut attributes = vec![
                Attribute {
                    tag: "ID".into(),
                    values: vec![id],
                },
                Attribute {
                    tag: "Name".into(),
                    values: vec![feature.label().into()],
                },
            ];
            if kind != feature.kind {
                attributes.push(Attribute {
                    tag: "cg2_original_type".into(),
                    values: vec![feature.kind.clone()],
                });
                warnings.push(format!("Feature {} type {:?} is exported as generic sequence_feature; original type remains in cg2_original_type, without guessing ontology equivalence.",index+1,feature.kind));
            }
            for key in ["gene", "product", "note"] {
                let values = feature
                    .qualifiers
                    .iter()
                    .filter(|(k, _)| k == key)
                    .filter_map(|(_, v)| v.clone())
                    .collect::<Vec<_>>();
                if !values.is_empty() {
                    attributes.push(Attribute {
                        tag: if key == "note" {
                            "Note".into()
                        } else {
                            format!("cg2_{key}")
                        },
                        values,
                    });
                }
            }
            if kind == "CDS" {
                let values = feature
                    .qualifiers
                    .iter()
                    .filter(|(k, _)| k == "codon_start")
                    .collect::<Vec<_>>();
                if values.len() != 1 {
                    return Err(unsupported(
                        "generic CDS GFF3 export requires one explicit codon_start; phase is never guessed",
                    ));
                }
                let first = values[0]
                    .1
                    .as_ref()
                    .ok_or_else(|| bad("codon_start lacks a value"))?
                    .parse::<u8>()
                    .ok()
                    .filter(|n| (1..=3).contains(n))
                    .ok_or_else(|| bad("invalid CDS codon_start"))?
                    - 1;
                if reverse {
                    parts.reverse();
                }
                let mut length = 0;
                for (i, &(a, b)) in parts.iter().enumerate() {
                    let phase = if i == 0 {
                        first
                    } else {
                        ((3 - (length % 3 + 3 - first as usize) % 3) % 3) as u8
                    };
                    let end = b;
                    push(
                        &mut out,
                        &gff_line(
                            seqid,
                            "Collagene",
                            kind,
                            a + 1,
                            end,
                            ".",
                            if reverse { "-" } else { "+" },
                            Some(phase),
                            &attributes,
                        ),
                    )?;
                    length += b - a;
                }
                if origin {
                    return Err(unsupported(
                        "origin-crossing generic CDS phase export requires retained source row evidence",
                    ));
                }
            } else {
                if origin {
                    let a = parts[0].0;
                    let end = reference.sequence.len() + parts[1].1;
                    push(
                        &mut out,
                        &gff_line(
                            seqid,
                            "Collagene",
                            kind,
                            a + 1,
                            end,
                            ".",
                            if reverse { "-" } else { "+" },
                            None,
                            &attributes,
                        ),
                    )?;
                } else {
                    for &(a, b) in &parts {
                        push(
                            &mut out,
                            &gff_line(
                                seqid,
                                "Collagene",
                                kind,
                                a + 1,
                                b,
                                ".",
                                if reverse { "-" } else { "+" },
                                None,
                                &attributes,
                            ),
                        )?;
                    }
                }
            }
            warnings.push(format!("Feature {} exported using current coordinates, label and explicit CDS phase where supplied; other qualifiers, source IDs/relationships and rendering colors are not inferred.",index+1));
        }
    }
    if format == AnnotationFormat::Gff3 && !reference.features.is_empty() {
        // Reparse our bounded output to catch missing Parent IDs, duplicated IDs,
        // illegal points/phases and contradictory row grouping before returning it.
        let _ = parse_gff(reference, &out, seqid)?;
    }
    if format == AnnotationFormat::Gtf && !reference.features.is_empty() {
        let _ = gtf::parse(reference, &out, seqid)?;
    }
    warnings.sort();
    warnings.dedup();
    Ok(ExportReport {
        text: out,
        warnings,
    })
}

/// Retrieve exact retained annotation source independently of subsequent edits.
/// Its hash is integrity evidence; no source coordinates are asserted as current.
pub fn original_annotation_source(
    reference: &Molecule,
    source_blake3: &str,
) -> Result<Vec<u8>, Error> {
    reference.validate()?;
    if source_blake3.len() != 64
        || !source_blake3
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
    {
        return Err(bad("invalid retained annotation source hash"));
    }
    let key = format!("annotation_import:{source_blake3}:original_base64");
    let encoded = reference
        .metadata
        .get(&key)
        .ok_or_else(|| bad("annotation source is not retained in this document"))?;
    if encoded.len() > MAX_ANNOTATION_BYTES.div_ceil(3) * 4 {
        return Err(Error::TooLarge);
    }
    let bytes = STANDARD
        .decode(encoded)
        .map_err(|_| bad("invalid retained annotation source encoding"))?;
    if bytes.len() > MAX_ANNOTATION_BYTES || blake3::hash(&bytes).to_hex().as_str() != source_blake3
    {
        return Err(bad("retained annotation source hash mismatch"));
    }
    Ok(bytes)
}
