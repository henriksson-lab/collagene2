//! Frame external GenBank records before the rich header/qualifier reader runs.
//! Feature tables are parsed independently and omitted from the private reader
//! copy. CONTIG gets a safe placeholder; original bytes remain unchanged.
use super::{MAX_IMPORT_BYTES, MAX_LOCATION_WORK, MAX_QUALIFIER_LINES, feature_qualifiers};
use collagene_model::{Error, Feature, MAX_SEQUENCE_LEN};

const MAX_RECORDS: usize = 1000;
const MAX_FEATURES: usize = 100_000;

fn invalid(message: &str) -> Error {
    Error::Format(format!("GenBank location framing: {message}"))
}

pub(super) struct ReaderInput<'a> {
    text: String,
    records: Vec<Record<'a>>,
}

struct Record<'a> {
    original: &'a str,
    features: Vec<FeatureLocation>,
}

struct FeatureLocation {
    kind: String,
    expression: String,
    qualifiers: Vec<(String, Option<String>)>,
}

struct RawRecord<'a> {
    original: &'a str,
    features: Vec<RawFeature>,
}

struct RawFeature {
    kind: String,
    expression: String,
    qualifier_lines: Vec<String>,
}

struct Pending {
    kind: Option<String>,
    expression: String,
    indent: usize,
}

#[derive(Default)]
struct Budget {
    features: usize,
    location_work: usize,
    qualifiers: usize,
    qualifier_lines: usize,
}

impl Budget {
    fn qualifier(&mut self, feature: &mut RawFeature, line: &str) -> Result<(), Error> {
        self.qualifier_lines = self
            .qualifier_lines
            .checked_add(1)
            .filter(|n| *n <= MAX_QUALIFIER_LINES)
            .ok_or(Error::TooLarge)?;
        feature.qualifier_lines.push(line.to_owned());
        Ok(())
    }
    fn location(&mut self, text: &str) -> Result<(), Error> {
        collagene_model::annotations::insdc::check_text(text)?;
        collagene_model::annotations::insdc::validate_for_genbank(text)?;
        // Conservative structural work bound across every feature and CONTIG,
        // in addition to the parser's per-expression byte/node/depth limits.
        let work = 1 + text
            .bytes()
            .filter(|b| matches!(b, b',' | b'(' | b':'))
            .count();
        self.location_work = self
            .location_work
            .checked_add(work)
            .filter(|n| *n <= MAX_LOCATION_WORK)
            .ok_or(Error::TooLarge)?;
        Ok(())
    }
}

impl Pending {
    fn append(&mut self, value: &str) -> Result<(), Error> {
        let value = value.trim_matches(' ');
        if self
            .expression
            .len()
            .checked_add(value.len())
            .filter(|n| *n <= 65_536)
            .is_none()
        {
            return Err(Error::TooLarge);
        }
        self.expression.push_str(value);
        Ok(())
    }

    fn finish(self, record: &mut RawRecord<'_>, budget: &mut Budget) -> Result<(), Error> {
        budget.location(&self.expression)?;
        if let Some(kind) = self.kind {
            budget.features = budget
                .features
                .checked_add(1)
                .filter(|n| *n <= MAX_FEATURES)
                .ok_or(Error::TooLarge)?;
            record.features.push(RawFeature {
                kind,
                expression: self.expression,
                qualifier_lines: Vec::new(),
            });
        }
        Ok(())
    }
}

impl<'a> ReaderInput<'a> {
    pub(super) fn new(text: &'a str) -> Result<Self, Error> {
        if text.len() > MAX_IMPORT_BYTES {
            return Err(Error::TooLarge);
        }
        let mut prepared = Self {
            text: String::with_capacity(text.len()),
            records: Vec::new(),
        };
        let mut budget = Budget::default();
        let mut start = 0;
        let mut offset = 0;
        for line in text.split_inclusive('\n') {
            offset += line.len();
            if line.starts_with("//") && line.trim_end() == "//" {
                prepared.record(&text[start..offset], &mut budget)?;
                start = offset;
            }
        }
        if start < text.len() && !text[start..].trim().is_empty() {
            prepared.record(&text[start..], &mut budget)?;
        }
        Ok(prepared)
    }

    fn record(&mut self, original: &'a str, budget: &mut Budget) -> Result<(), Error> {
        if self.records.len() >= MAX_RECORDS {
            return Err(Error::TooLarge);
        }
        let mut record = RawRecord {
            original,
            features: Vec::new(),
        };
        let mut in_features = false;
        let mut feature_indent = None;
        let mut pending: Option<Pending> = None;
        for raw in original.split_inclusive('\n') {
            let body = raw.strip_suffix('\n').unwrap_or(raw);
            let body = body.strip_suffix('\r').unwrap_or(body);
            let indent = body.bytes().take_while(|b| *b == b' ').count();
            if body.trim_matches(' ').is_empty() {
                if pending.is_none() && in_features {
                    if let Some(feature) = record.features.last_mut() {
                        budget.qualifier(feature, "")?;
                    }
                } else if pending.is_none() {
                    self.text.push_str(raw);
                }
                continue;
            }
            if indent == 0 {
                if let Some(location) = pending.take() {
                    location.finish(&mut record, budget)?;
                }
                in_features = body.starts_with("FEATURES")
                    && matches!(body.as_bytes().get(8), Some(b' ' | b'\t'));
                feature_indent = None;
                if let Some(value) = body.strip_prefix("CONTIG").filter(|s| s.starts_with(' ')) {
                    let prefix = body.len() - value.trim_start_matches(' ').len();
                    let mut location = Pending {
                        kind: None,
                        expression: String::new(),
                        indent: prefix,
                    };
                    location.append(&body[prefix..])?;
                    pending = Some(location);
                    self.placeholder(&body[..prefix], raw);
                    continue;
                }
                if body.starts_with("LOCUS") && matches!(body.as_bytes().get(5), Some(b' ' | b'\t'))
                {
                    // The upstream LOCUS parser reads u64 then casts to usize.
                    // Reject an oversized declared count before that cast on WASM.
                    if let Some(length) = body.split_ascii_whitespace().nth(2) {
                        if let Ok(length) = length.parse::<u64>() {
                            if length > MAX_SEQUENCE_LEN as u64 {
                                return Err(Error::TooLarge);
                            }
                        }
                    }
                }
                self.text.push_str(raw);
                continue;
            }
            if let Some(location) = &mut pending {
                if indent >= location.indent {
                    let value = &body[location.indent..];
                    if location.kind.is_some() && value.trim_start().starts_with('/') {
                        pending
                            .take()
                            .expect("pending location")
                            .finish(&mut record, budget)?;
                        budget.qualifier(
                            record
                                .features
                                .last_mut()
                                .ok_or_else(|| invalid("qualifier has no feature"))?,
                            value,
                        )?;
                    } else {
                        location.append(value)?;
                    }
                    continue;
                }
                pending
                    .take()
                    .expect("pending location")
                    .finish(&mut record, budget)?;
            }
            if in_features && !body.trim().is_empty() && feature_indent.is_none_or(|n| indent < n) {
                let value = &body[indent..];
                let key_end = value
                    .bytes()
                    .position(|b| b == b' ')
                    .ok_or_else(|| invalid("feature location is missing"))?;
                let kind = &value[..key_end];
                if !super::feature_kind_valid(kind) {
                    return Err(invalid(
                        "invalid feature key or qualifier without a feature",
                    ));
                }
                let prefix =
                    indent + key_end + value[key_end..].bytes().take_while(|b| *b == b' ').count();
                feature_indent = Some(prefix);
                let mut location = Pending {
                    kind: Some(kind.into()),
                    expression: String::new(),
                    indent: prefix,
                };
                location.append(&body[prefix..])?;
                pending = Some(location);
            } else if in_features {
                let column = feature_indent.ok_or_else(|| invalid("qualifier has no feature"))?;
                budget.qualifier(
                    record
                        .features
                        .last_mut()
                        .ok_or_else(|| invalid("qualifier has no feature"))?,
                    body.get(column..)
                        .ok_or_else(|| invalid("invalid qualifier indentation"))?,
                )?;
            } else {
                self.text.push_str(raw);
            }
        }
        if let Some(location) = pending {
            location.finish(&mut record, budget)?;
        }
        let features = record
            .features
            .into_iter()
            .map(|raw| {
                let qualifiers = feature_qualifiers::parse(&raw.qualifier_lines)?;
                budget.qualifiers = budget
                    .qualifiers
                    .checked_add(qualifiers.len())
                    .filter(|n| *n <= 1_000_000)
                    .ok_or(Error::TooLarge)?;
                Ok(FeatureLocation {
                    kind: raw.kind,
                    expression: raw.expression,
                    qualifiers,
                })
            })
            .collect::<Result<Vec<_>, Error>>()?;
        self.records.push(Record {
            original: record.original,
            features,
        });
        Ok(())
    }

    fn placeholder(&mut self, prefix: &str, raw: &str) {
        self.text.push_str(prefix);
        self.text.push('1');
        self.text
            .push_str(if raw.ends_with("\r\n") { "\r\n" } else { "\n" });
    }

    pub(super) fn text(&self) -> &str {
        &self.text
    }

    pub(super) fn original(&self, index: usize) -> Result<&'a str, Error> {
        self.records
            .get(index)
            .map(|r| r.original)
            .ok_or_else(|| invalid("original record is missing"))
    }

    pub(super) fn features(&self, index: usize, length: usize) -> Result<Vec<Feature>, Error> {
        let record = self
            .records
            .get(index)
            .ok_or_else(|| invalid("feature record is missing"))?;
        record
            .features
            .iter()
            .map(|feature| {
                // Reparse with actual length: unresolved mixed remote expressions
                // can still contain local coordinates that must remain bounded.
                let location = collagene_model::annotations::insdc::parse_location(
                    &feature.expression,
                    length,
                )?;
                let needs_review = feature.qualifiers.iter().any(|(key, value)| {
                    key == "collagene_review_required" && value.as_deref() == Some("true")
                });
                Ok(Feature {
                    kind: feature.kind.clone(),
                    location,
                    needs_review,
                    qualifiers: feature
                        .qualifiers
                        .iter()
                        .map(|(key, value)| super::import_qualifier(key, value))
                        .collect(),
                })
            })
            .collect()
    }

    pub(super) fn record_count(&self) -> usize {
        self.records.len()
    }
}
