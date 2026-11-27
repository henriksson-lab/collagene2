//! Independently implemented sequence screening, not PTM/structure prediction.
//! Supplied CDR intervals are documentary context, never validated numbering.
use collagene_model::{Alphabet, Molecule};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
pub const ENGINE: &str = "protein-liabilities-1";
pub const MAX_RESIDUES: usize = 1_000_000;
pub const MAX_HITS: usize = 100_000;
pub const MAX_CDRS: usize = 256;
const MAX_SOURCE_BYTES: usize = 32 << 20;
const MAX_OPTIONS_BYTES: usize = 512 << 10;
const AA: &[u8] = b"ACDEFGHIKLMNPQRSTVWYOU";
const LU: &str = "https://pmc.ncbi.nlm.nih.gov/articles/PMC6343770/";
const LAP: &str = "https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011881";
const BENCHLING: &str = "https://help.benchling.com/hc/en-us/articles/9684246077325-Annotate-antibody-specific-information-on-protein-sequences";
pub type Result<T> = std::result::Result<T, Error>;
#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub enum Error {
    #[error("invalid liability-screen input: {0}")]
    Invalid(&'static str),
    #[error("liability-screen resource limit: {0}")]
    Limit(&'static str),
    #[error("liability-screen source/options/report binding changed")]
    Binding,
    #[error("invalid source model: {0}")]
    Source(String),
    #[error("liability screening cancelled")]
    Cancelled,
}
macro_rules! data {($($t:tt)*)=>{#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)] #[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))] #[serde(deny_unknown_fields)] $($t)*};}
data! {pub enum CdrKind {Cdr1,Cdr2,Cdr3}}
data! {pub struct CdrRegion {pub domain_id:String,pub kind:CdrKind,pub start:usize,pub end:usize}}
data! {pub struct ProvidedCdrs {
    pub source_revision:String,
    /// User-supplied description only, not a certified IMGT/Kabat/Chothia mapping.
    pub definition_label:String,pub provenance:String,pub regions:Vec<CdrRegion>,
}}
data! {pub struct Options {pub max_hits:usize,pub provided_cdrs:Option<ProvidedCdrs>}}
data! {pub enum RuleId {CanonicalAsnDeamidation,CanonicalAspIsomerization,NGlycosylationSequon,MethionineOxidation,TryptophanOxidation,ProvidedCdrAsparagineReview}}
data! {pub enum Certainty {CertainPattern,PossiblePattern}}
data! {pub enum CdrMembership {Unassessed,OutsideProvided,Inside{region_index:usize}}}
data! {pub struct Hit {pub rule:RuleId,pub start:usize,pub end:usize,pub residues:String,pub certainty:Certainty,pub cdr:CdrMembership}}
data! {pub struct Counts {pub rule:RuleId,pub certain:usize,pub possible:usize}}
data! {pub struct Report {
    pub engine:String,pub catalog_revision:String,pub source_revision:String,pub context_revision:String,
    pub source_length:usize,pub options:Options,pub counts:Vec<Counts>,pub hits:Vec<Hit>,
    pub hits_truncated:bool,pub ambiguous_residues:usize,pub stop_symbols:usize,
    /// There is deliberately no automatic domain/CDR or numbering computation here.
    pub cdr_rule_unassessed:bool,
}}
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct Rule {
    pub id: RuleId,
    pub pattern: &'static str,
    pub masks: Vec<u32>,
    pub provided_cdr_only: bool,
    pub reference: &'static str,
}
fn mask(b: u8) -> u32 {
    match b {
        b'B' => mask(b'D') | mask(b'N'),
        b'Z' => mask(b'E') | mask(b'Q'),
        b'J' => mask(b'I') | mask(b'L'),
        b'X' => (1 << AA.len()) - 1,
        b'*' => 0,
        _ => AA.iter().position(|&a| a == b).map(|i| 1 << i).unwrap_or(0),
    }
}
fn group(b: &[u8]) -> u32 {
    b.iter().fold(0, |m, &b| m | mask(b))
}
pub fn rules() -> Vec<Rule> {
    use RuleId::*;
    vec![
        Rule {
            id: CanonicalAsnDeamidation,
            pattern: "N[GSNTH]",
            masks: vec![mask(b'N'), group(b"GSNTH")],
            provided_cdr_only: false,
            reference: LU,
        },
        Rule {
            id: CanonicalAspIsomerization,
            pattern: "D[GSDTH]",
            masks: vec![mask(b'D'), group(b"GSDTH")],
            provided_cdr_only: false,
            reference: LU,
        },
        Rule {
            id: NGlycosylationSequon,
            pattern: "N[^P][ST]",
            masks: vec![mask(b'N'), mask(b'X') & !mask(b'P'), group(b"ST")],
            provided_cdr_only: false,
            reference: LAP,
        },
        Rule {
            id: MethionineOxidation,
            pattern: "M",
            masks: vec![mask(b'M')],
            provided_cdr_only: false,
            reference: LAP,
        },
        Rule {
            id: TryptophanOxidation,
            pattern: "W",
            masks: vec![mask(b'W')],
            provided_cdr_only: false,
            reference: LAP,
        },
        Rule {
            id: ProvidedCdrAsparagineReview,
            pattern: "N within supplied CDR",
            masks: vec![mask(b'N')],
            provided_cdr_only: true,
            reference: BENCHLING,
        },
    ]
}
fn hash<T: Serialize>(domain: &[u8], value: &T, limit: usize) -> Result<String> {
    struct Sink {
        h: blake3::Hasher,
        n: usize,
        limit: usize,
    }
    impl std::io::Write for Sink {
        fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
            self.n = self
                .n
                .checked_add(b.len())
                .filter(|&n| n <= self.limit)
                .ok_or_else(|| std::io::Error::other("screen JSON bound"))?;
            self.h.update(b);
            Ok(b.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut out = Sink {
        h: blake3::Hasher::new(),
        n: 0,
        limit,
    };
    out.h.update(domain);
    serde_json::to_writer(&mut out, value)
        .map_err(|_| Error::Limit("serialized source/options"))?;
    Ok(out.h.finalize().to_hex().to_string())
}
pub fn catalog_revision() -> String {
    hash(b"protein-liability-catalog-1\0",&(ENGINE,rules(),AA,mask(b'B'),mask(b'Z'),mask(b'J'),mask(b'X'),"anchor-in-provided-CDR; adjacent source context may extend outside region; stop blocks windows; no probability/chemistry/numbering authority"),65536).expect("fixed rule catalog")
}
pub fn source_revision(m: &Molecule) -> Result<String> {
    if m.sequence.alphabet() != Alphabet::Protein {
        return Err(Error::Invalid("explicit protein source required"));
    }
    if m.sequence.len() > MAX_RESIDUES {
        return Err(Error::Limit("source residues"));
    }
    collagene_model::bounds::ModelBudget::default()
        .molecule(m)
        .map_err(|e| Error::Source(e.to_string()))?;
    m.validate().map_err(|e| Error::Source(e.to_string()))?;
    hash(b"", m, MAX_SOURCE_BYTES)
}
fn text(s: &str, max: usize) -> Result<()> {
    if s.trim().is_empty() || s.len() > max || s.chars().any(char::is_control) {
        Err(Error::Invalid(
            "empty, oversized or control-containing label",
        ))
    } else {
        Ok(())
    }
}
fn options_context<'a>(
    revision: &str,
    length: usize,
    o: &'a Options,
) -> Result<Vec<(usize, &'a CdrRegion)>> {
    if o.max_hits > MAX_HITS {
        return Err(Error::Limit("retained hits"));
    }
    let mut ordered = vec![];
    if let Some(c) = &o.provided_cdrs {
        if c.source_revision != revision {
            return Err(Error::Binding);
        }
        text(&c.definition_label, 128)?;
        text(&c.provenance, 1024)?;
        if c.regions.len() > MAX_CDRS {
            return Err(Error::Limit("provided CDR regions"));
        }
        let mut unique = BTreeSet::new();
        for (i, r) in c.regions.iter().enumerate() {
            text(&r.domain_id, 64)?;
            let kind = match r.kind {
                CdrKind::Cdr1 => 1,
                CdrKind::Cdr2 => 2,
                CdrKind::Cdr3 => 3,
            };
            if r.start >= r.end || r.end > length || !unique.insert((&r.domain_id, kind)) {
                return Err(Error::Invalid("invalid or repeated provided CDR"));
            }
            ordered.push((i, r));
        }
        ordered.sort_by_key(|(_, r)| r.start);
        if ordered.windows(2).any(|p| p[0].1.end > p[1].1.start) {
            return Err(Error::Invalid("provided CDRs overlap"));
        }
    }
    hash(b"protein-liability-options-1\0", o, MAX_OPTIONS_BYTES)?;
    Ok(ordered)
}
fn prepare<'a>(m: &Molecule, o: &'a Options) -> Result<(String, Vec<(usize, &'a CdrRegion)>)> {
    let revision = source_revision(m)?;
    let ordered = options_context(&revision, m.sequence.len(), o)?;
    if ordered
        .iter()
        .any(|(_, r)| m.sequence.as_str().as_bytes()[r.start..r.end].contains(&b'*'))
    {
        return Err(Error::Invalid("provided CDR crosses a stop"));
    }
    Ok((revision, ordered))
}
fn membership(o: &Options, ordered: &[(usize, &CdrRegion)], position: usize) -> CdrMembership {
    let i = ordered.partition_point(|(_, r)| r.end <= position);
    if o.provided_cdrs.is_none() {
        CdrMembership::Unassessed
    } else if let Some((index, r)) = ordered.get(i).filter(|(_, r)| r.start <= position) {
        let _ = r;
        CdrMembership::Inside {
            region_index: *index,
        }
    } else {
        CdrMembership::OutsideProvided
    }
}
impl RuleId {
    pub fn label(&self) -> &'static str {
        match self {
            Self::CanonicalAsnDeamidation => "Asn canonical-motif review",
            Self::CanonicalAspIsomerization => "Asp canonical-motif review",
            Self::NGlycosylationSequon => "N-glycosylation sequon",
            Self::MethionineOxidation => "Methionine oxidation review",
            Self::TryptophanOxidation => "Tryptophan oxidation review",
            Self::ProvidedCdrAsparagineReview => "Asparagine in supplied CDR",
        }
    }
}
pub fn request_revision(m: &Molecule, o: &Options) -> Result<String> {
    let (source, _) = prepare(m, o)?;
    hash(
        b"protein-liability-request-1\0",
        &(ENGINE, source, catalog_revision(), o),
        MAX_OPTIONS_BYTES + 65536,
    )
}
fn classify(bytes: &[u8], required: &[u32]) -> Option<Certainty> {
    let mut certain = true;
    for (&b, &required) in bytes.iter().zip(required) {
        let candidates = mask(b);
        if candidates & required == 0 {
            return None;
        }
        certain &= candidates & !required == 0;
    }
    Some(if certain {
        Certainty::CertainPattern
    } else {
        Certainty::PossiblePattern
    })
}
pub fn screen(m: &Molecule, o: &Options, cancel: &mut impl FnMut() -> bool) -> Result<Report> {
    let (source_revision, ordered) = prepare(m, o)?;
    let catalog_revision = catalog_revision();
    let context_revision = hash(
        b"protein-liability-request-1\0",
        &(ENGINE, &source_revision, &catalog_revision, o),
        MAX_OPTIONS_BYTES + 65536,
    )?;
    let rules = rules();
    let bytes = m.sequence.as_str().as_bytes();
    let mut counts = rules
        .iter()
        .map(|r| Counts {
            rule: r.id.clone(),
            certain: 0,
            possible: 0,
        })
        .collect::<Vec<_>>();
    let mut hits = Vec::with_capacity(o.max_hits.min(4096));
    let mut total = 0;
    let mut ambiguous_residues = 0;
    let mut stop_symbols = 0;
    for i in 0..bytes.len() {
        if i % 1024 == 0 && cancel() {
            return Err(Error::Cancelled);
        }
        ambiguous_residues += usize::from(matches!(bytes[i], b'B' | b'J' | b'X' | b'Z'));
        stop_symbols += usize::from(bytes[i] == b'*');
        let membership = membership(o, &ordered, i);
        for (k, r) in rules.iter().enumerate() {
            if r.provided_cdr_only && !matches!(membership, CdrMembership::Inside { .. }) {
                continue;
            }
            let end = i + r.masks.len();
            if end > bytes.len() {
                continue;
            }
            let Some(certainty) = classify(&bytes[i..end], &r.masks) else {
                continue;
            };
            match certainty {
                Certainty::CertainPattern => counts[k].certain += 1,
                Certainty::PossiblePattern => counts[k].possible += 1,
            };
            total += 1;
            if hits.len() < o.max_hits {
                hits.push(Hit {
                    rule: r.id.clone(),
                    start: i,
                    end,
                    residues: m.sequence.as_str()[i..end].into(),
                    certainty,
                    cdr: membership.clone(),
                });
            }
        }
    }
    if cancel() {
        return Err(Error::Cancelled);
    }
    Ok(Report {
        engine: ENGINE.into(),
        catalog_revision,
        source_revision,
        context_revision,
        source_length: bytes.len(),
        options: o.clone(),
        counts,
        hits_truncated: total > hits.len(),
        hits,
        ambiguous_residues,
        stop_symbols,
        cdr_rule_unassessed: o
            .provided_cdrs
            .as_ref()
            .is_none_or(|c| c.regions.is_empty()),
    })
}
/// Validate one independently supplied match against the actual immutable source.
/// This establishes a literal candidate pattern only, never a measured PTM.
pub fn validate_hit(m: &Molecule, o: &Options, h: &Hit) -> Result<()> {
    let (_, ordered) = prepare(m, o)?;
    let rules = rules();
    let r = rules
        .iter()
        .find(|r| r.id == h.rule)
        .ok_or(Error::Binding)?;
    let cdr = membership(o, &ordered, h.start);
    let end = h.start.checked_add(r.masks.len()).ok_or(Error::Binding)?;
    let bytes = m
        .sequence
        .as_str()
        .as_bytes()
        .get(h.start..end)
        .ok_or(Error::Binding)?;
    if end != h.end
        || bytes != h.residues.as_bytes()
        || classify(bytes, &r.masks) != Some(h.certainty.clone())
        || cdr != h.cdr
        || (r.provided_cdr_only && !matches!(cdr, CdrMembership::Inside { .. }))
    {
        return Err(Error::Binding);
    }
    Ok(())
}
impl Report {
    /// Source-free transport consistency only. Saved evidence must still use
    /// validate_for against independently retained original source/options.
    pub fn validate_binding(&self, expected: &str) -> Result<()> {
        let source = &self.source_revision;
        let rules = rules();
        if self.engine != ENGINE
            || self.catalog_revision != catalog_revision()
            || source.len() != 64
            || !source
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
            || self.source_length == 0
            || self.source_length > MAX_RESIDUES
            || self.counts.len() != rules.len()
            || self.hits.len() > MAX_HITS
            || self.ambiguous_residues > self.source_length
            || self.stop_symbols > self.source_length - self.ambiguous_residues
            || self.context_revision != expected
            || hash(
                b"protein-liability-request-1\0",
                &(ENGINE, source, &self.catalog_revision, &self.options),
                MAX_OPTIONS_BYTES + 65536,
            )? != self.context_revision
        {
            return Err(Error::Binding);
        }
        let ordered = options_context(source, self.source_length, &self.options)?;
        let unassessed = self
            .options
            .provided_cdrs
            .as_ref()
            .is_none_or(|c| c.regions.is_empty());
        if self.cdr_rule_unassessed != unassessed {
            return Err(Error::Binding);
        }
        let mut total = 0;
        for (c, r) in self.counts.iter().zip(&rules) {
            if c.rule != r.id
                || c.certain > self.source_length
                || c.possible > self.source_length - c.certain
            {
                return Err(Error::Binding);
            }
            let count = c.certain + c.possible;
            if r.provided_cdr_only
                && count > ordered.iter().map(|(_, r)| r.end - r.start).sum::<usize>()
            {
                return Err(Error::Binding);
            }
            total += count;
        }
        if self.hits.len() != total.min(self.options.max_hits)
            || self.hits_truncated != (total > self.options.max_hits)
        {
            return Err(Error::Binding);
        }
        let mut previous = None;
        let mut observed = BTreeMap::new();
        let mut retained = vec![(0usize, 0usize); rules.len()];
        for h in &self.hits {
            let k = rules
                .iter()
                .position(|r| r.id == h.rule)
                .ok_or(Error::Binding)?;
            let r = &rules[k];
            let end = h.start.checked_add(r.masks.len()).ok_or(Error::Binding)?;
            if h.end != end
                || end > self.source_length
                || h.residues.len() != r.masks.len()
                || !h.residues.is_ascii()
                || previous.is_some_and(|p| p >= (h.start, k))
                || h.cdr != membership(&self.options, &ordered, h.start)
                || (r.provided_cdr_only && !matches!(h.cdr, CdrMembership::Inside { .. }))
                || classify(h.residues.as_bytes(), &r.masks) != Some(h.certainty.clone())
            {
                return Err(Error::Binding);
            }
            previous = Some((h.start, k));
            for (i, b) in h.residues.bytes().enumerate() {
                if observed.insert(h.start + i, b).is_some_and(|old| old != b) {
                    return Err(Error::Binding);
                }
            }
            match h.certainty {
                Certainty::CertainPattern => retained[k].0 += 1,
                Certainty::PossiblePattern => retained[k].1 += 1,
            };
        }
        if observed.len() > self.source_length - self.stop_symbols
            || observed
                .values()
                .filter(|&&b| matches!(b, b'B' | b'Z' | b'J' | b'X'))
                .count()
                > self.ambiguous_residues
        {
            return Err(Error::Binding);
        }
        for ((certain, possible), c) in retained.iter().zip(&self.counts) {
            if *certain > c.certain
                || *possible > c.possible
                || (!self.hits_truncated && (*certain != c.certain || *possible != c.possible))
            {
                return Err(Error::Binding);
            }
        }
        Ok(())
    }
    /// All matches, counts, context and truncation are freshly reconstructed.
    pub fn validate_for(&self, m: &Molecule, o: &Options) -> Result<()> {
        if self.engine != ENGINE
            || self.hits.len() > MAX_HITS
            || self.counts.len() != rules().len()
            || screen(m, o, &mut || false)? != *self
        {
            return Err(Error::Binding);
        }
        Ok(())
    }
}
