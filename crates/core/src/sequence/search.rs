//! Bounded, source-bound searches. Coordinates are zero-based half-open.
use crate::sequence::analysis::{base_mask, starter_enzymes};
use collagene_model::*;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum SearchMode {
    All,
    Sequence,
    Annotations,
    Primers,
    Enzymes,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct SearchOptions {
    pub query: String,
    pub mode: SearchMode,
    pub both_strands: bool,
    pub case_sensitive: bool,
    pub whole_word: bool,
    pub enzymes: Vec<String>,
    pub max_hits: usize,
    pub max_rows: usize,
    pub max_work: usize,
}
impl Default for SearchOptions {
    fn default() -> Self {
        Self {
            query: String::new(),
            mode: SearchMode::All,
            both_strands: true,
            case_sensitive: false,
            whole_word: false,
            enzymes: vec![],
            max_hits: 1000,
            max_rows: 100_000,
            max_work: 10_000_000,
        }
    }
}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum SearchKind {
    Sequence,
    Annotation,
    Primer,
    Enzyme,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct SearchSpan {
    pub start: usize,
    pub end: usize,
    pub reverse: bool,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct SearchHit {
    pub id: String,
    pub kind: SearchKind,
    pub label: String,
    pub parts: Vec<SearchSpan>,
    pub feature_index: Option<usize>,
    pub enzyme: Option<String>,
    pub top_cut: Option<usize>,
    pub bottom_cut: Option<usize>,
    pub certain: bool,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct SearchReport {
    pub source_length: usize,
    pub source_feature_count: usize,
    pub source_revision: String,
    pub engine_version: String,
    pub options: SearchOptions,
    pub hits: Vec<SearchHit>,
    pub complete: bool,
    pub work: usize,
    pub rows: usize,
    pub warnings: Vec<String>,
}
impl SearchOptions {
    pub fn validate(&self) -> Result<(), Error> {
        if self.query.trim().is_empty() || self.query.len() > 4096 {
            return Err(Error::Format(
                "search query must contain 1..4096 bytes".into(),
            ));
        }
        if self.max_hits == 0
            || self.max_hits > 10_000
            || self.max_rows == 0
            || self.max_rows > 100_000
            || self.max_work == 0
            || self.max_work > 50_000_000
            || self.enzymes.len() > 100
            || self.enzymes.iter().any(|s| s.len() > 64)
        {
            return Err(Error::Unsupported(
                "search limits exceed supported bounds".into(),
            ));
        }
        Ok(())
    }
}
impl SearchReport {
    pub fn validate(&self) -> Result<(), Error> {
        self.options.validate()?;
        let bad = || Error::Format("invalid search report bounds or identity".into());
        if self.source_length == 0
            || self.source_length > collagene_model::MAX_SEQUENCE_LEN
            || self.source_feature_count > 100_000
            || self.source_revision.len() != 64
            || !self.source_revision.bytes().all(|b| b.is_ascii_hexdigit())
            || self.engine_version.is_empty()
            || self.engine_version.len() > 256
            || self.hits.len() > self.options.max_hits
            || self.work > self.options.max_work
            || self.rows > self.options.max_rows
            || self.warnings.len() > 100
        {
            return Err(bad());
        }
        let mut ids = std::collections::BTreeSet::new();
        for h in &self.hits {
            if h.label.len() > 4096
                || h.parts.len() > 4096
                || h.parts
                    .iter()
                    .any(|p| p.start >= p.end || p.end > self.source_length)
                || h.feature_index
                    .is_some_and(|i| i >= self.source_feature_count)
                || h.top_cut.is_some_and(|i| i > self.source_length)
                || h.bottom_cut.is_some_and(|i| i > self.source_length)
                || !ids.insert(&h.id)
            {
                return Err(bad());
            }
            if matches!(h.kind, SearchKind::Annotation | SearchKind::Primer)
                != h.feature_index.is_some()
                || (h.kind == SearchKind::Enzyme) != h.enzyme.is_some()
                || (h.kind != SearchKind::Enzyme && (h.top_cut.is_some() || h.bottom_cut.is_some()))
                || (matches!(h.kind, SearchKind::Sequence | SearchKind::Enzyme)
                    && h.parts.is_empty())
            {
                return Err(bad());
            }
            let mut original = h.clone();
            original.id.clear();
            let expected = blake3::hash(
                &serde_json::to_vec(&(
                    &self.source_revision,
                    &self.engine_version,
                    &self.options,
                    &original,
                ))
                .map_err(|_| bad())?,
            )
            .to_hex()
            .to_string();
            if h.id != expected {
                return Err(bad());
            }
        }
        Ok(())
    }
}

fn text_match(text: &str, o: &SearchOptions, q: &str) -> bool {
    if text.is_empty() {
        return false;
    }
    let s = if o.case_sensitive {
        text.to_owned()
    } else {
        text.to_lowercase()
    };
    if s.len() < q.len() {
        return false;
    }
    s.match_indices(q).any(|(i, _)| {
        !o.whole_word
            || (!s[..i]
                .chars()
                .next_back()
                .is_some_and(|c| c.is_alphanumeric() || c == '_')
                && !s[i + q.len()..]
                    .chars()
                    .next()
                    .is_some_and(|c| c.is_alphanumeric() || c == '_'))
    })
}
fn location_parts(
    l: &Location,
    reverse: bool,
    out: &mut Vec<SearchSpan>,
    work: &mut usize,
    limit: usize,
) -> bool {
    if *work == limit || out.len() == 4096 {
        return false;
    }
    *work += 1;
    match l {
        Location::Range { start, end, .. } => out.push(SearchSpan {
            start: *start,
            end: *end,
            reverse,
        }),
        Location::Complement(inner) => return location_parts(inner, !reverse, out, work, limit),
        Location::Join(parts) | Location::Order(parts) => {
            if reverse {
                for p in parts.iter().rev() {
                    if !location_parts(p, reverse, out, work, limit) {
                        return false;
                    }
                }
            } else {
                for p in parts {
                    if !location_parts(p, reverse, out, work, limit) {
                        return false;
                    }
                }
            }
        }
        Location::Unresolved(_) => {}
    }
    true
}
fn exact_location(l: &Location) -> bool {
    match l {
        Location::Range { before, after, .. } => !*before && !*after,
        Location::Complement(l) => exact_location(l),
        Location::Join(parts) => parts.iter().all(exact_location),
        _ => false,
    }
}

fn spans(start: usize, len: usize, n: usize, reverse: bool) -> Vec<SearchSpan> {
    let mut parts = vec![SearchSpan {
        start,
        end: (start + len).min(n),
        reverse,
    }];
    if start + len > n {
        parts.push(SearchSpan {
            start: 0,
            end: start + len - n,
            reverse,
        });
    }
    if reverse {
        parts.reverse();
    }
    parts
}
fn push(r: &mut SearchReport, mut h: SearchHit) -> bool {
    if r.hits.len() == r.options.max_hits {
        r.complete = false;
        return false;
    }
    h.id = blake3::hash(
        &serde_json::to_vec(&(&r.source_revision, &r.engine_version, &r.options, &h))
            .expect("serializable search hit"),
    )
    .to_hex()
    .to_string();
    r.hits.push(h);
    true
}
fn hit(kind: SearchKind, label: String, parts: Vec<SearchSpan>) -> SearchHit {
    SearchHit {
        id: String::new(),
        kind,
        label,
        parts,
        feature_index: None,
        enzyme: None,
        top_cut: None,
        bottom_cut: None,
        certain: true,
    }
}
fn sequence_scan(
    r: &mut SearchReport,
    m: &Molecule,
    pattern: &Sequence,
    enzyme: Option<&crate::sequence::analysis::Enzyme>,
) -> Result<(), Error> {
    let n = m.sequence.len();
    let k = pattern.len();
    if k > n {
        return Ok(());
    }
    let mut patterns = vec![(pattern.clone(), false)];
    if r.options.both_strands && m.sequence.alphabet() != Alphabet::Protein {
        let rc = pattern.reverse_complement()?;
        if rc != *pattern {
            patterns.push((rc, true));
        }
    }
    let limit = if m.topology == Topology::Circular {
        n
    } else {
        n - k + 1
    };
    for (p, reverse) in patterns {
        for start in 0..limit {
            let mut certain = true;
            let mut matched = true;
            for (j, b) in p.as_str().bytes().enumerate() {
                if r.work == r.options.max_work {
                    r.complete = false;
                    return Ok(());
                }
                r.work += 1;
                let a = m.sequence.as_str().as_bytes()[(start + j) % n];
                if m.sequence.alphabet() == Alphabet::Protein {
                    if a != b {
                        matched = false;
                        break;
                    }
                } else {
                    let (a, b) = (base_mask(a), base_mask(b));
                    if a & b == 0 {
                        matched = false;
                        break;
                    }
                    certain &= a & !b == 0;
                }
            }
            if matched {
                let mut h = hit(
                    if enzyme.is_some() {
                        SearchKind::Enzyme
                    } else {
                        SearchKind::Sequence
                    },
                    enzyme.map_or_else(|| r.options.query.clone(), |e| e.name.clone()),
                    spans(start, k, n, reverse),
                );
                h.certain = certain;
                if let Some(e) = enzyme {
                    h.enzyme = Some(e.name.clone());
                    let (t, b) = if reverse {
                        (
                            start as i64 + k as i64 - e.bottom_cut as i64,
                            start as i64 + k as i64 - e.top_cut as i64,
                        )
                    } else {
                        (
                            start as i64 + e.top_cut as i64,
                            start as i64 + e.bottom_cut as i64,
                        )
                    };
                    let pos = |x: i64| {
                        if m.topology == Topology::Circular {
                            Some(x.rem_euclid(n as i64) as usize)
                        } else if (0..=n as i64).contains(&x) {
                            Some(x as usize)
                        } else {
                            None
                        }
                    };
                    h.top_cut = pos(t);
                    h.bottom_cut = pos(b);
                }
                if !push(r, h) {
                    return Ok(());
                }
            }
        }
    }
    Ok(())
}
pub fn search_document(m: &Molecule, o: &SearchOptions) -> Result<SearchReport, Error> {
    m.validate()?;
    o.validate()?;
    let normalized_query = if o.case_sensitive {
        o.query.clone()
    } else {
        o.query.to_lowercase()
    };
    let mut r = SearchReport {
        source_length: m.sequence.len(),
        source_feature_count: m.features.len(),
        source_revision: m.revision(),
        engine_version: crate::verification::ENGINE_VERSION.into(),
        options: o.clone(),
        hits: vec![],
        complete: true,
        work: 0,
        rows: 0,
        warnings: vec![],
    };
    if matches!(o.mode, SearchMode::All | SearchMode::Sequence) {
        if m.sequence.alphabet() == Alphabet::Protein && o.both_strands {
            r.warnings.push(
                "Protein sequence search is forward-only; reverse complement is undefined.".into(),
            );
        }
        match Sequence::new(m.sequence.alphabet(), &o.query) {
            Ok(p) => sequence_scan(&mut r, m, &p, None)?,
            Err(e) => {
                if o.mode == SearchMode::Sequence {
                    return Err(e);
                } else {
                    r.warnings.push("Query is not a sequence in this document alphabet; sequence source skipped.".into())
                }
            }
        }
    }
    if matches!(
        o.mode,
        SearchMode::All | SearchMode::Annotations | SearchMode::Primers
    ) {
        'features: for (i, f) in m.features.iter().enumerate() {
            if r.rows == o.max_rows || r.work == o.max_work {
                r.complete = false;
                break;
            }
            r.rows += 1;
            r.work += 1;
            let primer = f.kind.eq_ignore_ascii_case("primer_bind");
            if o.mode == SearchMode::Primers && !primer {
                continue;
            }
            let mut remaining = o.max_work - r.work;
            if f.kind.len() > remaining {
                r.work = o.max_work;
                r.complete = false;
                break;
            }
            remaining -= f.kind.len();
            for (key, value) in &f.qualifiers {
                let cost = 1usize
                    .saturating_add(key.len())
                    .saturating_add(value.as_ref().map_or(0, String::len));
                if cost > remaining {
                    r.work = o.max_work;
                    r.complete = false;
                    break 'features;
                }
                remaining -= cost;
            }
            r.work = o.max_work - remaining;
            if text_match(&f.kind, o, &normalized_query)
                || f.qualifiers.iter().any(|(k, v)| {
                    text_match(k, o, &normalized_query)
                        || v.as_deref()
                            .is_some_and(|s| text_match(s, o, &normalized_query))
                })
            {
                let mut parts = vec![];
                if !location_parts(&f.location, false, &mut parts, &mut r.work, o.max_work) {
                    r.complete = false;
                    if parts.len() == 4096 {
                        let warning="An annotation exceeded the 4096-part coordinate limit and was omitted.".to_string();
                        if !r.warnings.contains(&warning) {
                            r.warnings.push(warning);
                        }
                    }

                    continue;
                }
                let mut h = hit(
                    if primer {
                        SearchKind::Primer
                    } else {
                        SearchKind::Annotation
                    },
                    f.label().chars().take(512).collect(),
                    parts,
                );
                h.feature_index = Some(i);
                h.certain = !f.needs_review && exact_location(&f.location);
                if !push(&mut r, h) {
                    break;
                }
            }
        }
    }
    if matches!(o.mode, SearchMode::All | SearchMode::Enzymes) {
        if m.sequence.alphabet() != Alphabet::Dna {
            if o.mode == SearchMode::Enzymes {
                return Err(Error::Unsupported("enzyme search requires DNA".into()));
            }
            r.warnings
                .push("Enzyme search requires DNA; enzyme source skipped.".into());
        } else {
            let catalog = starter_enzymes();
            for name in &o.enzymes {
                if !catalog.iter().any(|e| e.name.eq_ignore_ascii_case(name)) {
                    return Err(Error::Unsupported(format!("unknown enzyme {name}")));
                }
            }
            for e in catalog {
                if (!o.enzymes.is_empty()
                    && !o.enzymes.iter().any(|n| n.eq_ignore_ascii_case(&e.name)))
                    || !(text_match(&e.name, o, &normalized_query)
                        || text_match(&e.recognition, o, &normalized_query))
                {
                    continue;
                }
                let p = Sequence::dna(&e.recognition)?;
                sequence_scan(&mut r, m, &p, Some(&e))?;
            }
        }
    }
    if !r.complete {
        r.warnings.push("Search stopped at a configured hit, row or work bound; results are incomplete across the requested sources.".into());
    }
    Ok(r)
}
