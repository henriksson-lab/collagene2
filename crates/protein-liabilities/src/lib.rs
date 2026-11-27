//! Independently implemented sequence screening, not PTM/structure prediction.
//! Supplied CDR intervals are documentary context, never validated numbering.
use collagene_model::{Alphabet, Molecule};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
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
fn prepare<'a>(m: &Molecule, o: &'a Options) -> Result<(String, Vec<(usize, &'a CdrRegion)>)> {
    let revision = source_revision(m)?;
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
            if r.start >= r.end
                || r.end > m.sequence.len()
                || !unique.insert((&r.domain_id, kind))
                || m.sequence.as_str().as_bytes()[r.start..r.end].contains(&b'*')
            {
                return Err(Error::Invalid(
                    "invalid, repeated or stop-crossing provided CDR",
                ));
            }
            ordered.push((i, r));
        }
        ordered.sort_by_key(|(_, r)| r.start);
        if ordered.windows(2).any(|p| p[0].1.end > p[1].1.start) {
            return Err(Error::Invalid("provided CDRs overlap"));
        }
    }
    hash(b"protein-liability-options-1\0", o, MAX_OPTIONS_BYTES)?;
    Ok((revision, ordered))
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
    let mut cursor = 0;
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
        while cursor < ordered.len() && ordered[cursor].1.end <= i {
            cursor += 1;
        }
        let membership = if o.provided_cdrs.is_none() {
            CdrMembership::Unassessed
        } else if cursor < ordered.len() && ordered[cursor].1.start <= i {
            CdrMembership::Inside {
                region_index: ordered[cursor].0,
            }
        } else {
            CdrMembership::OutsideProvided
        };
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
impl Report {
    /// All matches, counts, context and truncation are freshly reconstructed.
    pub fn validate_for(&self, m: &Molecule, o: &Options) -> Result<()> {
        if self.engine != ENGINE
            || self.hits.len() > MAX_HITS
            || self.counts.len() != 6
            || screen(m, o, &mut || false)? != *self
        {
            return Err(Error::Binding);
        }
        Ok(())
    }
}
