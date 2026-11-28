//! Independent GTF row exchange, deliberately not a transcript assembler.
use super::*;

fn tag(s: &str) -> Result<(), Error> {
    if s.is_empty()
        || s.len() > 128
        || !s.bytes().enumerate().all(|(i, b)| {
            b.is_ascii_alphabetic()
                || b == b'_'
                || i > 0 && (b.is_ascii_digit() || b".:-".contains(&b))
        })
    {
        return Err(bad("invalid GTF attribute tag"));
    }
    Ok(())
}
fn token(s: &str, label: &str) -> Result<(), Error> {
    text_token(s, label)?;
    if s.chars().any(char::is_whitespace) {
        return Err(bad(format!("GTF {label} cannot contain whitespace")));
    }
    Ok(())
}
fn attributes(text: &str) -> Result<Vec<Attribute>, Error> {
    if text == "." {
        return Ok(vec![]);
    }
    if text.trim().is_empty() {
        return Err(bad(
            "missing GTF attribute field; use '.' for an explicit missing value",
        ));
    }
    let bytes = text.as_bytes();
    let mut at = 0;
    let mut result = vec![];
    while at < bytes.len() {
        while bytes.get(at) == Some(&b' ') {
            at += 1;
        }
        if at == bytes.len() {
            break;
        }
        if result.len() >= 1000 {
            return Err(Error::TooLarge);
        }
        let start = at;
        while at < bytes.len() && bytes[at] != b' ' {
            at += 1;
        }
        let name = &text[start..at];
        tag(name)?;
        if at == bytes.len() {
            return Err(bad("GTF attribute lacks a separated value"));
        }
        while bytes.get(at) == Some(&b' ') {
            at += 1;
        }
        let value = if bytes.get(at) == Some(&b'"') {
            at += 1;
            let mut value = Vec::new();
            let mut closed = false;
            while at < bytes.len() {
                let b = bytes[at];
                at += 1;
                if b == b'"' {
                    closed = true;
                    break;
                }
                if b == b'\\' {
                    let escaped = *bytes.get(at).ok_or_else(|| bad("truncated GTF escape"))?;
                    at += 1;
                    if !matches!(escaped, b'"' | b'\\') {
                        return Err(unsupported(
                            "GTF supports only explicit quote/backslash string escapes; other escapes require review",
                        ));
                    }
                    value.push(escaped);
                } else {
                    value.push(b);
                }
                if value.len() > 8192 {
                    return Err(Error::TooLarge);
                }
            }
            if !closed {
                return Err(bad("unterminated GTF attribute quotation"));
            }
            String::from_utf8(value).map_err(|_| bad("GTF attribute is not UTF-8"))?
        } else {
            // Common GTF2 numeric attributes are accepted but normalized to
            // quoted text. GFF3 tag=value and unquoted text are not GTF strings.
            let start = at;
            while at < bytes.len() && bytes[at] != b';' && bytes[at] != b' ' {
                at += 1;
            }
            let value = &text[start..at];
            if value.is_empty()
                || value.len() > 8192
                || value.parse::<f64>().ok().is_none_or(|v| !v.is_finite())
            {
                return Err(bad(
                    "GTF string values must be quoted; only finite numeric values may be unquoted",
                ));
            }
            value.to_owned()
        };
        if value.chars().any(char::is_control) {
            return Err(bad("GTF attribute contains a control character"));
        }
        while bytes.get(at) == Some(&b' ') {
            at += 1;
        }
        if bytes.get(at) != Some(&b';') {
            return Err(bad("GTF attributes require terminating semicolons"));
        }
        at += 1;
        if bytes.get(at).is_some_and(|b| *b != b' ') {
            return Err(bad("GTF attributes must be separated by spaces"));
        }
        result.push(Attribute {
            tag: name.into(),
            values: vec![value],
        });
    }
    Ok(result)
}
fn validate(row: &GffRow, reference: Option<&Molecule>) -> Result<(), Error> {
    seqid_token(&row.seqid)?;
    token(&row.seqid, "sequence ID")?;
    token(&row.source, "source")?;
    token(&row.kind, "feature type")?;
    if row.start == 0 || row.end < row.start {
        return Err(bad("GTF coordinates require 1 <= start <= end"));
    }
    if reference.is_some_and(|r| row.end > r.sequence.len()) {
        return Err(bad(
            "GTF interval exceeds the explicitly mapped reference; origin wrapping is not inferred",
        ));
    }
    if row.score != "." && row.score.parse::<f64>().ok().is_none_or(|v| !v.is_finite()) {
        return Err(bad("GTF score must be finite or missing"));
    }
    if !matches!(row.strand.as_str(), "+" | "-" | ".") || row.phase.is_some_and(|p| p > 2) {
        return Err(bad("invalid GTF strand/frame"));
    }
    if row.attributes.len() > 1000 {
        return Err(Error::TooLarge);
    }
    for a in &row.attributes {
        tag(&a.tag)?;
        if a.values.len() != 1
            || a.values[0].len() > 8192
            || a.values[0].chars().any(char::is_control)
        {
            return Err(bad("invalid retained GTF attribute"));
        }
    }
    Ok(())
}
pub(super) fn parse(
    reference: &Molecule,
    text: &str,
    seqid: &str,
) -> Result<(Vec<AnnotationProposal>, Vec<String>, Vec<String>), Error> {
    token(seqid, "sequence ID")?;
    let mut proposals = vec![];
    let mut skipped = BTreeSet::new();
    let mut directives = vec![];
    let mut rows = 0;
    for (index, line) in text.lines().enumerate() {
        let line = line.strip_suffix('\r').unwrap_or(line);
        if line.trim().is_empty() {
            continue;
        }
        if line.starts_with('#') || line.starts_with("track ") || line.starts_with("browser ") {
            if directives.len() >= MAX_DIRECTIVES {
                return Err(Error::TooLarge);
            }
            if line == "##FASTA" || line.starts_with("##gff-version 3") {
                return Err(unsupported(
                    "GFF3/embedded FASTA is not this GTF row profile",
                ));
            }
            directives.push(line.into());
            continue;
        }
        rows += 1;
        if rows > MAX_ANNOTATION_ROWS {
            return Err(Error::TooLarge);
        }
        let f = line.split('\t').collect::<Vec<_>>();
        if f.len() != 9 {
            return Err(bad(format!(
                "GTF line {} requires nine tab-separated fields",
                index + 1
            )));
        }
        let phase = match f[7] {
            "." => None,
            "0" => Some(0),
            "1" => Some(1),
            "2" => Some(2),
            _ => return Err(bad("GTF frame must be 0, 1, 2 or missing")),
        };
        let row = GffRow {
            line: index + 1,
            seqid: f[0].into(),
            source: f[1].into(),
            kind: f[2].into(),
            start: number(f[3])?,
            end: number(f[4])?,
            score: f[5].into(),
            strand: f[6].into(),
            phase,
            attributes: attributes(f[8])?,
        };
        validate(&row, None)?;
        if row.seqid == seqid {
            proposals.push(proposal(row, reference)?);
        } else {
            skipped.insert(row.seqid);
        }
    }
    Ok((proposals, skipped.into_iter().collect(), directives))
}
pub(super) fn proposal(row: GffRow, reference: &Molecule) -> Result<AnnotationProposal, Error> {
    validate(&row, Some(reference))?;
    let mut warnings=vec!["GTF rows are retained independently; gene_id/transcript_id attributes are recorded relationships, not reconstructed transcripts, splice validity or sequence verification.".into()];
    let range = Location::range(row.start - 1, row.end);
    let location = match row.strand.as_str() {
        "-" => Location::Complement(Box::new(range)),
        "+" => range,
        _ => {
            warnings.push("GTF strand is unspecified; coordinates remain visible but biological extraction is disabled.".into());
            Location::Order(vec![range])
        }
    };
    let label = [
        "cg2_label",
        "gene_name",
        "transcript_name",
        "transcript_id",
        "gene_id",
    ]
    .into_iter()
    .find_map(|k| scalar(&row, k))
    .unwrap_or_else(|| row.kind.clone());
    let mut qualifiers = vec![
        ("label".into(), Some(label)),
        ("gtf:source".into(), Some(row.source.clone())),
        ("gtf:score".into(), Some(row.score.clone())),
        ("gtf:strand".into(), Some(row.strand.clone())),
        (
            "gtf:phase".into(),
            Some(
                row.phase
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| ".".into()),
            ),
        ),
    ];
    for a in &row.attributes {
        qualifiers.push((
            format!("gtf:attribute:{}", a.tag),
            Some(a.values[0].clone()),
        ));
    }
    for key in ["gene_id", "transcript_id"] {
        if key == "transcript_id" && row.kind == "gene" {
            continue;
        }
        let values = row
            .attributes
            .iter()
            .filter(|a| a.tag == key)
            .collect::<Vec<_>>();
        if values.len() != 1 || values[0].values[0].is_empty() {
            warnings.push(format!("GTF {key} is missing, empty or repeated; no unique gene/transcript assignment is inferred."));
        }
    }
    if row.kind == "CDS" {
        if let Some(phase) = row.phase.filter(|_| row.strand != ".") {
            qualifiers.push(("codon_start".into(), Some((phase + 1).to_string())));
        } else {
            warnings.push("CDS phase or strand is unspecified; no current codon_start or translation is asserted.".into());
        }
    } else if row.phase.is_some() {
        warnings.push(
            "A non-CDS frame is preserved as source evidence, without inferring a coding sequence."
                .into(),
        );
    }
    let evidence = SourceRecord {
        format: AnnotationFormat::Gtf,
        kind: row.kind.clone(),
        location: location.clone(),
        gff: vec![row.clone()],
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
            kind: row.kind,
            location,
            qualifiers,
            needs_review: true,
        },
        source_rows: vec![row.line],
        warnings,
    })
}
fn quote(value: &str) -> String {
    format!("\"{}\"", value.replace('\\', "\\\\").replace('"', "\\\""))
}
pub(super) fn line(seqid: &str, row: &GffRow) -> Result<String, Error> {
    token(seqid, "sequence ID")?;
    validate(row, None)?;
    let mut size = seqid
        .len()
        .checked_add(row.source.len())
        .and_then(|n| n.checked_add(row.kind.len()))
        .and_then(|n| n.checked_add(row.score.len()))
        .and_then(|n| n.checked_add(128))
        .ok_or(Error::TooLarge)?;
    for a in &row.attributes {
        size = size
            .checked_add(a.tag.len())
            .and_then(|n| n.checked_add(5))
            .and_then(|n| n.checked_add(a.values[0].len()))
            .and_then(|n| {
                n.checked_add(
                    a.values[0]
                        .bytes()
                        .filter(|b| matches!(*b, b'"' | b'\\'))
                        .count(),
                )
            })
            .filter(|n| *n <= 64 * 1024)
            .ok_or(Error::TooLarge)?;
    }
    let attrs = if row.attributes.is_empty() {
        ".".into()
    } else {
        row.attributes
            .iter()
            .map(|a| format!("{} {};", a.tag, quote(&a.values[0])))
            .collect::<Vec<_>>()
            .join(" ")
    };
    Ok(format!(
        "{seqid}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{attrs}\n",
        row.source,
        row.kind,
        row.start,
        row.end,
        row.score,
        row.strand,
        row.phase
            .map(|p| p.to_string())
            .unwrap_or_else(|| ".".into())
    ))
}
/// Generic export requires explicit relationship IDs, never invented hierarchy.
pub(super) fn generic(
    feature: &Feature,
    reference: &Molecule,
    seqid: &str,
    index: usize,
    out: &mut String,
    warnings: &mut Vec<String>,
) -> Result<(), Error> {
    let (mut parts, reverse) = exact_parts(&feature.location)?;
    if parts.is_empty() || parts.windows(2).any(|w| w[0].1 > w[1].0) {
        return Err(unsupported(
            "GTF generic export requires reference-ordered nonoverlapping intervals; circular origin joins require separate reviewed modeling",
        ));
    }
    let mut attrs = vec![];
    for (key, value) in &feature.qualifiers {
        if let Some(key) = key.strip_prefix("gtf:attribute:") {
            tag(key)?;
            let value = value
                .as_ref()
                .ok_or_else(|| bad("GTF attribute lacks a value"))?;
            if attrs.len() >= 999 || value.len() > 8192 || value.chars().any(char::is_control) {
                return Err(Error::TooLarge);
            }
            attrs.push(Attribute {
                tag: key.into(),
                values: vec![value.clone()],
            });
        }
    }
    for key in ["gene_id", "transcript_id"] {
        if key == "transcript_id" && feature.kind == "gene" {
            continue;
        }
        let values = attrs.iter().filter(|a| a.tag == key).collect::<Vec<_>>();
        if values.len() != 1 || values[0].values[0].is_empty() {
            return Err(unsupported(format!(
                "GTF current-feature export requires one explicit gtf:attribute:{key}; synthetic gene/transcript identity is not invented"
            )));
        }
    }
    if feature.label().len() > 8192 || feature.label().chars().any(char::is_control) {
        return Err(Error::TooLarge);
    }
    attrs.push(Attribute {
        tag: "cg2_label".into(),
        values: vec![feature.label().into()],
    });
    let initial = if feature.kind == "CDS" {
        let values = feature
            .qualifiers
            .iter()
            .filter(|(k, _)| k == "codon_start")
            .collect::<Vec<_>>();
        if values.len() != 1 {
            return Err(unsupported(
                "GTF current CDS export requires one explicit codon_start; phase is not guessed",
            ));
        }
        Some(
            values[0]
                .1
                .as_deref()
                .and_then(|s| s.parse::<u8>().ok())
                .filter(|v| (1..=3).contains(v))
                .ok_or_else(|| bad("invalid CDS codon_start"))?
                - 1,
        )
    } else {
        None
    };
    if reverse {
        parts.reverse();
    }
    let mut length = 0;
    for &(a, b) in &parts {
        let phase = initial.map(|first| {
            if length == 0 {
                first
            } else {
                ((3 - (length % 3 + 3 - first as usize) % 3) % 3) as u8
            }
        });
        let row = GffRow {
            line: 0,
            seqid: seqid.into(),
            source: "Collagene".into(),
            kind: feature.kind.clone(),
            start: a + 1,
            end: b,
            score: ".".into(),
            strand: if reverse { "-" } else { "+" }.into(),
            phase,
            attributes: attrs.clone(),
        };
        validate(&row, Some(reference))?;
        push(out, &line(seqid, &row)?)?;
        length += b - a;
    }
    warnings.push(format!("Feature {} exports current intervals with explicit gene/transcript IDs and CDS phase where supplied. Other qualifiers and source score/phase/relationships are not inferred; native evidence remains separate. Joined parts remain rows, not a newly validated transcript model.",index+1));
    Ok(())
}
