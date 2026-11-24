//! Independently implemented single nucleotide record routing, from official ENA/DDBJ docs.
//! ENA: https://ena-docs.readthedocs.io/en/latest/retrieval/programmatic-access/browser-api.html
//! DDBJ: https://www.ddbj.nig.ac.jp/services/getentry-e.html
use super::*;

// Restrict this adapter to nucleotide record accessions, excluding project/sample/run,
// assembly accessions and WGS/TSA master records. No numeric UID or collection lookup.
pub(super) fn sequence_accession(value: &str) -> Option<String> {
    let id = accession(value, false)?;
    let base = id.split('.').next()?;
    let letters = base.bytes().take_while(u8::is_ascii_alphabetic).count();
    let digits = &base[letters..];
    if !digits.bytes().all(|b| b.is_ascii_digit())
        || !matches!(
            (letters, digits.len()),
            (1, 5) | (2, 6) | (2, 8) | (4, 8) | (4, 10) | (6, 9) | (6, 10)
        )
        || (letters >= 4 && digits.ends_with("000000"))
    {
        return None;
    }
    Some(id)
}

pub(super) fn resolve(url: &Url) -> Option<Result<Source, SourceError>> {
    let host = url.host_str()?;
    let ena = matches!(host, "www.ebi.ac.uk" | "ebi.ac.uk") && url.path().starts_with("/ena/");
    let ddbj = host == "getentry.ddbj.nig.ac.jp";
    if !ena && !ddbj {
        return None;
    }
    Some((|| {
        // Empty internal segments and encoded selectors are never normalized away.
        let path: Vec<_> = url
            .path()
            .trim_end_matches('/')
            .split('/')
            .skip(1)
            .collect();
        let q: Vec<_> = url.query_pairs().collect();
        let mut seen = std::collections::BTreeSet::new();
        if q.iter().any(|(k, _)| !seen.insert(k.as_ref())) {
            return Err(SourceError::UnsupportedSelection);
        }
        let id = if ena {
            if q.iter()
                .any(|(k, v)| k != "download" || !matches!(v.as_ref(), "true" | "false"))
            {
                return Err(SourceError::UnsupportedSelection);
            }
            match path.as_slice() {
                ["ena", "browser", "view", id] => *id,
                ["ena", "browser", "api", "embl" | "fasta", id] => *id,
                _ => return Err(SourceError::UnsupportedSelection),
            }
            .to_owned()
        } else {
            for (key, value) in &q {
                let allowed = match key.as_ref() {
                    "database" => matches!(value.as_ref(), "na" | "ddbj"),
                    "accession_number" => true,
                    "format" => matches!(value.as_ref(), "flatfile" | "fasta"),
                    "filetype" => matches!(value.as_ref(), "text" | "html"),
                    "limit" => value == "1",
                    _ => false,
                };
                if !allowed {
                    return Err(SourceError::UnsupportedSelection);
                }
            }
            match path.as_slice() {
                ["getentry", "na" | "ddbj", id]
                    if !q
                        .iter()
                        .any(|(k, _)| matches!(k.as_ref(), "database" | "accession_number")) =>
                {
                    (*id).to_owned()
                }
                ["getentry"] => q
                    .iter()
                    .find(|(k, _)| k == "accession_number")
                    .map(|(_, v)| v.to_string())
                    .ok_or(SourceError::UnsupportedSelection)?,
                _ => return Err(SourceError::UnsupportedSelection),
            }
        };
        let accession = sequence_accession(&id).ok_or(SourceError::UnsupportedSelection)?;
        Ok(if ena {
            Source::Ena { accession }
        } else {
            Source::Ddbj { accession }
        })
    })())
}

pub(super) fn download(source: &Source) -> Result<Url, SourceError> {
    let (value, ena) = match source {
        Source::Ena { accession } => (accession, true),
        Source::Ddbj { accession } => (accession, false),
        _ => return Err(SourceError::InvalidSource),
    };
    if sequence_accession(value).as_ref() != Some(value) {
        return Err(SourceError::InvalidSource);
    }
    Url::parse(&if ena {
        format!("https://www.ebi.ac.uk/ena/browser/api/embl/{value}")
    } else {
        format!("https://getentry.ddbj.nig.ac.jp/getentry/na/{value}?format=flatfile&filetype=text")
    })
    .map_err(|_| SourceError::InvalidSource)
}

pub(super) fn actual_embl(m: &Molecule) -> Result<String, SourceError> {
    let raw = m
        .metadata
        .get("embl.original_record")
        .ok_or(SourceError::InvalidSequenceData)?;
    let id: Vec<_> = raw
        .lines()
        .find_map(|l| l.strip_prefix("ID   "))
        .ok_or(SourceError::MissingAccessionVersion)?
        .split(';')
        .map(str::trim)
        .collect();
    let primary = raw
        .lines()
        .filter_map(|l| l.strip_prefix("AC   "))
        .flat_map(|l| l.split(';'))
        .map(str::trim)
        .find(|s| !s.is_empty())
        .and_then(sequence_accession)
        .filter(|s| !s.contains('.'))
        .ok_or(SourceError::MissingAccessionVersion)?;
    let sv: Vec<_> = raw
        .lines()
        .filter_map(|l| l.strip_prefix("SV   "))
        .collect();
    if sv.len() > 1 {
        return Err(SourceError::IdentifierMismatch);
    }
    let actual = if id.len() == 7 {
        if id[0] != primary {
            return Err(SourceError::IdentifierMismatch);
        }
        let version = id[1]
            .strip_prefix("SV ")
            .ok_or(SourceError::MissingAccessionVersion)?;
        sequence_accession(&format!("{primary}.{version}"))
            .ok_or(SourceError::MissingAccessionVersion)?
    } else {
        // Legacy ID is a locus name. Its parser's fallback version 1 is not evidence.
        sv.first()
            .and_then(|s| sequence_accession(s.trim()))
            .filter(|s| s.contains('.'))
            .ok_or(SourceError::MissingAccessionVersion)?
    };
    if actual.split('.').next() != Some(primary.as_str())
        || sv.first().is_some_and(|s| s.trim() != actual)
    {
        return Err(SourceError::IdentifierMismatch);
    }
    Ok(actual)
}

pub(super) fn validate_genbank_identity(m: &Molecule, actual: &str) -> Result<(), SourceError> {
    let raw = m
        .original_genbank
        .as_ref()
        .ok_or(SourceError::InvalidSequenceData)?;
    let primary = raw
        .lines()
        .find(|l| l.starts_with("ACCESSION "))
        .and_then(|l| l.split_ascii_whitespace().nth(1))
        .and_then(sequence_accession)
        .ok_or(SourceError::MissingAccessionVersion)?;
    if actual.split('.').next() != Some(primary.as_str()) {
        return Err(SourceError::IdentifierMismatch);
    }
    Ok(())
}
