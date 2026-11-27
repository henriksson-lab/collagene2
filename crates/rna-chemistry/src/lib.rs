//! Declared covalent RNA chemistry; no measured-state, cloning or activity authority.
mod formula;
use collagene_model::{Alphabet, Molecule, Topology};
pub use formula::{Formula, MASS_MODEL, Mass};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
pub const ENGINE: &str = "rna-chemistry-2";
pub const DECLARATION_VERSION: u16 = 2;
pub const MAX_BASES: usize = 1_000_000;
pub const MAX_SITES: usize = 100_000;
pub const MAX_LINKS: usize = 100_000;
pub const MAX_DECLARATION_BYTES: usize = 8 << 20;
const MAX_SOURCE_BYTES: usize = 32 << 20;
pub type Result<T> = std::result::Result<T, Error>;
#[derive(Debug, thiserror::Error, PartialEq, Eq)]
pub enum Error {
    #[error("invalid RNA declaration: {0}")]
    Invalid(&'static str),
    #[error("RNA resource bound exceeded: {0}")]
    Limit(&'static str),
    #[error("RNA source/declaration/report binding changed")]
    Binding,
    #[error("invalid source model: {0}")]
    Source(String),
}
macro_rules! data {($($t:tt)*)=>{
    #[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]
    #[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]
    #[serde(deny_unknown_fields)] $($t)*
};}
mod annotations;
pub use annotations::{
    AnnotationCandidate, AnnotationProfile, FeatureReview, annotation_candidates,
    annotation_profiles,
};
data! {pub struct Provenance {pub statement:String,pub reference:Option<String>}}
data! {pub enum Coverage {
    /// Omitted sites/bonds may contain other chemistry; quantities stay unavailable.
    Partial,
    /// Explicit hypothesis: all unlisted residues/bonds use the canonical model.
    ClosedModel,
}}
data! {pub enum Nucleoside {Canonical,N6Methyladenosine,FiveMethylcytidine,FiveMethyluridine,Pseudouridine,Inosine,Unknown{label:String}}}
data! {pub enum Sugar {Ribose,TwoPrimeOMethyl,TwoPrimeDeoxy,TwoPrimeFluoro,Unknown{label:String}}}
data! {pub enum Stereo {Rp,Sp,Unspecified}}
data! {pub enum Linkage {Phosphodiester,Phosphorothioate{stereo:Stereo},Unknown{label:String}}}
data! {pub enum FivePrime {Hydroxyl,Monophosphate,Triphosphate,Unknown{label:String}}}
data! {pub enum ThreePrime {Hydroxyl,Monophosphate,TwoThreeCyclicPhosphate,Unknown{label:String}}}
data! {pub enum Ends {Linear{five:FivePrime,three:ThreePrime},Circular}}
data! {pub struct Site {pub index:usize,pub nucleoside:Nucleoside,pub sugar:Sugar,pub provenance:Provenance}}
data! {pub struct Link {/// Stored-strand zero-based upstream residue; last position is the circular seam only.
pub after:usize,pub linkage:Linkage,pub provenance:Provenance}}
data! {pub struct Declaration {
    pub version:u16,pub source_revision:String,pub coverage:Coverage,pub ends:Ends,
    pub sites:Vec<Site>,pub links:Vec<Link>,pub annotation_reviews:Vec<FeatureReview>,pub provenance:Provenance,
}}
data! {pub enum IssueKind {PartialCoverage,AmbiguousReference,UnknownNucleoside,UnknownSugar,UnknownLinkage,UnknownFivePrime,UnknownThreePrime,RetainedSourceEvidence}}
data! {pub struct Issue {pub kind:IssueKind,pub count:usize,pub positions:Vec<usize>,pub positions_truncated:bool}}
data! {pub enum RetainedEvidence {Feature{index:usize},Metadata{key:String}}}
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub engine: String,
    pub catalog_revision: String,
    pub source_revision: String,
    pub declaration_revision: String,
    pub context_revision: String,
    pub source_length: usize,
    pub topology: Topology,
    pub declaration: Declaration,
    pub formula: Option<Formula>,
    pub mass: Option<Mass>,
    /// This is a complete declared formula model, never a measured-state certificate.
    pub formula_complete: bool,
    /// Additionally requires specified PS stereochemistry; mass cannot distinguish it.
    pub identity_complete: bool,
    pub unspecified_stereocenters: Vec<usize>,
    pub retained_evidence: Vec<RetainedEvidence>,
    pub unreviewed_evidence: Vec<RetainedEvidence>,
    pub issues: Vec<Issue>,
}
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct Profile {
    pub nucleoside: Nucleoside,
    pub canonical_parent: char,
    pub nucleoside_formula: Formula,
    pub reference: &'static str,
}
pub fn profiles() -> Vec<Profile> {
    use Nucleoside::*;
    [
        (
            N6Methyladenosine,
            'A',
            [11, 15, 5, 4, 0, 0, 0],
            "https://www.ebi.ac.uk/chebi/CHEBI:21891",
        ),
        (
            FiveMethylcytidine,
            'C',
            [10, 15, 3, 5, 0, 0, 0],
            "https://www.ebi.ac.uk/chebi/CHEBI:20607",
        ),
        (
            FiveMethyluridine,
            'U',
            [10, 14, 2, 6, 0, 0, 0],
            "https://www.ebi.ac.uk/chebi/CHEBI:45996",
        ),
        (
            Pseudouridine,
            'U',
            [9, 12, 2, 6, 0, 0, 0],
            "https://www.ebi.ac.uk/chebi/CHEBI:17802",
        ),
        (
            Inosine,
            'A',
            [10, 12, 4, 5, 0, 0, 0],
            "https://www.ebi.ac.uk/chebi/CHEBI:17596",
        ),
    ]
    .into_iter()
    .map(
        |(nucleoside, canonical_parent, [c, h, n, o, p, s, f], reference)| Profile {
            nucleoside,
            canonical_parent,
            nucleoside_formula: Formula {
                c,
                h,
                n,
                o,
                p,
                s,
                f,
            },
            reference,
        },
    )
    .collect()
}
fn sugar_delta(s: &Sugar) -> Option<[i64; 7]> {
    match s {
        Sugar::Ribose => Some([0; 7]),
        Sugar::TwoPrimeOMethyl => Some([1, 2, 0, 0, 0, 0, 0]),
        Sugar::TwoPrimeDeoxy => Some([0, 0, 0, -1, 0, 0, 0]),
        Sugar::TwoPrimeFluoro => Some([0, -1, 0, -1, 0, 0, 1]),
        Sugar::Unknown { .. } => None,
    }
}
fn label(value: &str, max: usize) -> Result<()> {
    if value.trim().is_empty() || value.len() > max || value.chars().any(char::is_control) {
        Err(Error::Invalid(
            "empty, oversized or control-containing label",
        ))
    } else {
        Ok(())
    }
}
fn provenance(p: &Provenance) -> Result<()> {
    label(&p.statement, 1024)?;
    if let Some(s) = &p.reference {
        label(s, 4096)?;
    }
    Ok(())
}
fn hash<T: Serialize>(domain: &[u8], value: &T, limit: usize) -> Result<String> {
    struct Sink {
        hash: blake3::Hasher,
        n: usize,
        limit: usize,
    }
    impl std::io::Write for Sink {
        fn write(&mut self, bytes: &[u8]) -> std::io::Result<usize> {
            self.n = self
                .n
                .checked_add(bytes.len())
                .filter(|&n| n <= self.limit)
                .ok_or_else(|| std::io::Error::other("RNA serialized bound"))?;
            self.hash.update(bytes);
            Ok(bytes.len())
        }
        fn flush(&mut self) -> std::io::Result<()> {
            Ok(())
        }
    }
    let mut s = Sink {
        hash: blake3::Hasher::new(),
        n: 0,
        limit,
    };
    s.hash.update(domain);
    serde_json::to_writer(&mut s, value)
        .map_err(|_| Error::Limit("serialized source/declaration"))?;
    Ok(s.hash.finalize().to_hex().to_string())
}
pub fn catalog_revision() -> String {
    // Pin neutral nucleosides, atom masses and all covalent deltas, not just names.
    hash(
        b"collagene-rna-catalog-1\0",
        &(
            profiles(),
            annotation_profiles(),
            formula::NUCLEOSIDES,
            formula::MONO,
            formula::AVERAGE,
            formula::BOND,
            formula::PHOSPHATE,
            formula::TRIPHOSPHATE,
            sugar_delta(&Sugar::TwoPrimeOMethyl),
            sugar_delta(&Sugar::TwoPrimeDeoxy),
            sugar_delta(&Sugar::TwoPrimeFluoro),
            formula::PHOSPHOROTHIOATE,
            MASS_MODEL,
        ),
        65536,
    )
    .expect("fixed catalog bounds")
}
pub fn source_revision(source: &Molecule) -> Result<String> {
    if source.sequence.alphabet() != Alphabet::Rna {
        return Err(Error::Invalid("explicit RNA source required"));
    }
    if source.sequence.len() > MAX_BASES {
        return Err(Error::Limit("source bases"));
    }
    collagene_model::bounds::ModelBudget::default()
        .molecule(source)
        .map_err(|e| Error::Source(e.to_string()))?;
    source
        .validate()
        .map_err(|e| Error::Source(e.to_string()))?;
    // Same compact serializer as Molecule::revision, with a bounded streaming hash.
    hash(b"", source, MAX_SOURCE_BYTES)
}
pub fn declaration_revision(source: &Molecule, d: &Declaration) -> Result<String> {
    let revision = source_revision(source)?;
    if d.version != DECLARATION_VERSION {
        return Err(Error::Invalid("unsupported declaration version"));
    }
    if revision != d.source_revision {
        return Err(Error::Binding);
    }
    if d.sites.len() > MAX_SITES || d.links.len() > MAX_LINKS {
        return Err(Error::Limit("sites/links"));
    }
    annotations::validate_reviews(source, d)?;
    provenance(&d.provenance)?;
    if !matches!(
        (&d.ends, source.topology),
        (Ends::Linear { .. }, Topology::Linear) | (Ends::Circular, Topology::Circular)
    ) {
        return Err(Error::Invalid("terminal declaration contradicts topology"));
    }
    let mut sites = BTreeSet::new();
    let bases = source.sequence.as_str().as_bytes();
    let catalog = profiles();
    for s in &d.sites {
        if s.index >= bases.len() || !sites.insert(s.index) {
            return Err(Error::Invalid(
                "duplicate or out-of-range residue declaration",
            ));
        }
        provenance(&s.provenance)?;
        match &s.nucleoside {
            Nucleoside::Canonical => {}
            Nucleoside::Unknown { label: v } => label(v, 256)?,
            modified => {
                let p = catalog
                    .iter()
                    .find(|p| &p.nucleoside == modified)
                    .ok_or(Error::Invalid("unknown catalog profile"))?;
                if p.canonical_parent as u8 != bases[s.index] {
                    return Err(Error::Invalid(
                        "modification parent contradicts exact reference base",
                    ));
                }
            }
        }
        if let Sugar::Unknown { label: v } = &s.sugar {
            label(v, 256)?;
        }
    }
    let bonds = if source.topology == Topology::Circular {
        bases.len()
    } else {
        bases.len() - 1
    };
    let mut links = BTreeSet::new();
    for l in &d.links {
        if l.after >= bonds || !links.insert(l.after) {
            return Err(Error::Invalid(
                "duplicate, terminal or out-of-range linkage",
            ));
        }
        provenance(&l.provenance)?;
        if let Linkage::Unknown { label: v } = &l.linkage {
            label(v, 256)?;
        }
    }
    if let Ends::Linear { five, three } = &d.ends {
        if let FivePrime::Unknown { label: v } = five {
            label(v, 256)?;
        }
        if let ThreePrime::Unknown { label: v } = three {
            label(v, 256)?;
        }
        if matches!(three, ThreePrime::TwoThreeCyclicPhosphate) {
            if let Some(s) = d.sites.iter().find(|s| s.index == bases.len() - 1) {
                if matches!(
                    s.sugar,
                    Sugar::TwoPrimeOMethyl | Sugar::TwoPrimeDeoxy | Sugar::TwoPrimeFluoro
                ) {
                    return Err(Error::Invalid(
                        "2',3' cyclic phosphate requires an available 2' hydroxyl",
                    ));
                }
            }
        }
    }
    hash(b"collagene-rna-declaration-2\0", d, MAX_DECLARATION_BYTES)
}
pub fn request_revision(source: &Molecule, d: &Declaration) -> Result<String> {
    let declaration = declaration_revision(source, d)?;
    hash(
        b"collagene-rna-request-2\0",
        &(ENGINE, &d.source_revision, declaration, catalog_revision()),
        65536,
    )
}
fn retained_evidence(source: &Molecule) -> Vec<RetainedEvidence> {
    let mut out = vec![];
    for (index, f) in source.features.iter().enumerate() {
        if ["modified_base", "modified_nucleotide", "modified_residue"]
            .iter()
            .any(|k| f.kind.eq_ignore_ascii_case(k))
            || f.qualifiers
                .iter()
                .any(|(k, _)| k.eq_ignore_ascii_case("mod_base"))
        {
            out.push(RetainedEvidence::Feature { index });
        }
    }
    for key in source.metadata.keys() {
        if [
            "rna:modifications",
            "rna:modified",
            "rna_modifications",
            "modifications",
            "tracked_modifications",
        ]
        .iter()
        .any(|k| key.eq_ignore_ascii_case(k))
        {
            out.push(RetainedEvidence::Metadata { key: key.clone() });
        }
    }
    out
}
fn issue(out: &mut Vec<Issue>, kind: IssueKind, positions: impl Iterator<Item = usize>) {
    let mut preview = vec![];
    let mut count = 0;
    for p in positions {
        count += 1;
        if preview.len() < 100 {
            preview.push(p);
        }
    }
    if count > 0 {
        out.push(Issue {
            kind,
            count,
            positions_truncated: count > preview.len(),
            positions: preview,
        });
    }
}
pub fn analyze(source: &Molecule, d: &Declaration) -> Result<Report> {
    let declaration_revision = declaration_revision(source, d)?;
    let catalog_revision = catalog_revision();
    let context_revision = hash(
        b"collagene-rna-request-2\0",
        &(
            ENGINE,
            &d.source_revision,
            &declaration_revision,
            &catalog_revision,
        ),
        65536,
    )?;
    let bases = source.sequence.as_str().as_bytes();
    let mut issues = vec![];
    if d.coverage == Coverage::Partial {
        issues.push(Issue {
            kind: IssueKind::PartialCoverage,
            count: 1,
            positions: vec![],
            positions_truncated: false,
        });
    }
    issue(
        &mut issues,
        IssueKind::AmbiguousReference,
        bases
            .iter()
            .enumerate()
            .filter(|(_, b)| !matches!(**b, b'A' | b'C' | b'G' | b'U'))
            .map(|(i, _)| i),
    );
    issue(
        &mut issues,
        IssueKind::UnknownNucleoside,
        d.sites
            .iter()
            .filter(|s| matches!(s.nucleoside, Nucleoside::Unknown { .. }))
            .map(|s| s.index),
    );
    issue(
        &mut issues,
        IssueKind::UnknownSugar,
        d.sites
            .iter()
            .filter(|s| matches!(s.sugar, Sugar::Unknown { .. }))
            .map(|s| s.index),
    );
    issue(
        &mut issues,
        IssueKind::UnknownLinkage,
        d.links
            .iter()
            .filter(|s| matches!(s.linkage, Linkage::Unknown { .. }))
            .map(|s| s.after),
    );
    if let Ends::Linear { five, three } = &d.ends {
        if matches!(five, FivePrime::Unknown { .. }) {
            issues.push(Issue {
                kind: IssueKind::UnknownFivePrime,
                count: 1,
                positions: vec![0],
                positions_truncated: false,
            });
        }
        if matches!(three, ThreePrime::Unknown { .. }) {
            issues.push(Issue {
                kind: IssueKind::UnknownThreePrime,
                count: 1,
                positions: vec![bases.len() - 1],
                positions_truncated: false,
            });
        }
    }
    let retained_evidence = retained_evidence(source);
    let reviewed = d
        .annotation_reviews
        .iter()
        .map(|r| r.feature_index)
        .collect::<BTreeSet<_>>();
    let unreviewed_evidence = retained_evidence
        .iter()
        .filter(|e| match e {
            RetainedEvidence::Feature { index } => !reviewed.contains(index),
            RetainedEvidence::Metadata { .. } => true,
        })
        .cloned()
        .collect::<Vec<_>>();
    if !unreviewed_evidence.is_empty() {
        issues.push(Issue {
            kind: IssueKind::RetainedSourceEvidence,
            count: unreviewed_evidence.len(),
            positions: vec![],
            positions_truncated: false,
        });
    }
    let unspecified_stereocenters = d
        .links
        .iter()
        .filter(|l| {
            matches!(
                l.linkage,
                Linkage::Phosphorothioate {
                    stereo: Stereo::Unspecified
                }
            )
        })
        .map(|l| l.after)
        .collect::<Vec<_>>();
    let formula_complete = issues.is_empty();
    let identity_complete = formula_complete && unspecified_stereocenters.is_empty();
    let formula = if formula_complete {
        let mut sum = formula::Sum::default();
        for &base in bases {
            sum.add(
                formula::nucleoside(base).ok_or(Error::Invalid("unresolved reference base"))?,
                1,
            )?;
        }
        let bonds = if source.topology == Topology::Circular {
            bases.len()
        } else {
            bases.len() - 1
        };
        sum.add(formula::BOND, bonds)?;
        let profiles = profiles();
        for s in &d.sites {
            if s.nucleoside != Nucleoside::Canonical {
                let p = profiles
                    .iter()
                    .find(|p| p.nucleoside == s.nucleoside)
                    .ok_or(Error::Invalid("unresolved nucleoside"))?;
                let base = formula::nucleoside(bases[s.index])
                    .ok_or(Error::Invalid("unresolved reference"))?;
                let modified = p.nucleoside_formula.array();
                let delta = std::array::from_fn(|i| modified[i] as i64 - base[i]);
                sum.add(delta, 1)?;
            }
            sum.add(
                sugar_delta(&s.sugar).ok_or(Error::Invalid("unresolved sugar"))?,
                1,
            )?;
        }
        for l in &d.links {
            if matches!(l.linkage, Linkage::Phosphorothioate { .. }) {
                sum.add(formula::PHOSPHOROTHIOATE, 1)?;
            }
        }
        if let Ends::Linear { five, three } = &d.ends {
            match five {
                FivePrime::Hydroxyl => {}
                FivePrime::Monophosphate => sum.add(formula::PHOSPHATE, 1)?,
                FivePrime::Triphosphate => sum.add(formula::TRIPHOSPHATE, 1)?,
                FivePrime::Unknown { .. } => return Err(Error::Invalid("unresolved 5' end")),
            }
            match three {
                ThreePrime::Hydroxyl => {}
                ThreePrime::Monophosphate => sum.add(formula::PHOSPHATE, 1)?,
                ThreePrime::TwoThreeCyclicPhosphate => sum.add(formula::BOND, 1)?,
                ThreePrime::Unknown { .. } => return Err(Error::Invalid("unresolved 3' end")),
            }
        }
        Some(sum.finish()?)
    } else {
        None
    };
    let mass = formula.map(Mass::for_formula);
    Ok(Report {
        engine: ENGINE.into(),
        catalog_revision,
        source_revision: d.source_revision.clone(),
        declaration_revision,
        context_revision,
        source_length: bases.len(),
        topology: source.topology,
        declaration: d.clone(),
        formula,
        mass,
        formula_complete,
        identity_complete,
        unspecified_stereocenters,
        retained_evidence,
        unreviewed_evidence,
        issues,
    })
}
impl Report {
    /// Response binding and bounded shape consistency, without the original bases.
    /// This cannot establish chemical correctness: accepting a saved result still
    /// requires `validate_for` against its independently retained source.
    pub fn validate_binding(&self, expected: &str) -> Result<()> {
        let d = &self.declaration;
        if self.engine != ENGINE
            || d.version != DECLARATION_VERSION
            || self.catalog_revision != catalog_revision()
            || self.source_revision != d.source_revision
            || self.source_revision.len() != 64
            || !self
                .source_revision
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
            || self.source_length == 0
            || self.source_length > MAX_BASES
            || d.sites.len() > MAX_SITES
            || d.links.len() > MAX_LINKS
            || d.annotation_reviews.len() > 4096
            || self.retained_evidence.len() > 8192
            || self.unreviewed_evidence.len() > self.retained_evidence.len()
            || hash(b"collagene-rna-declaration-2\0", d, MAX_DECLARATION_BYTES)?
                != self.declaration_revision
            || hash(
                b"collagene-rna-request-2\0",
                &(
                    ENGINE,
                    &self.source_revision,
                    &self.declaration_revision,
                    &self.catalog_revision,
                ),
                65536,
            )? != self.context_revision
            || self.context_revision != expected
            || self.formula_complete != self.issues.is_empty()
            || self.formula_complete != self.formula.is_some()
            || self.identity_complete
                != (self.formula_complete && self.unspecified_stereocenters.is_empty())
            || self.mass != self.formula.map(Mass::for_formula)
        {
            return Err(Error::Binding);
        }
        let mut sites = BTreeSet::new();
        let mut links = BTreeSet::new();
        let mut reviews = BTreeSet::new();
        let bonds = if self.topology == Topology::Circular {
            self.source_length
        } else {
            self.source_length - 1
        };
        if d.sites
            .iter()
            .any(|s| s.index >= self.source_length || !sites.insert(s.index))
            || d.links
                .iter()
                .any(|l| l.after >= bonds || !links.insert(l.after))
            || d.annotation_reviews
                .iter()
                .any(|r| !reviews.insert(r.feature_index) || r.site_entry >= d.sites.len())
            || !matches!(
                (&d.ends, self.topology),
                (Ends::Linear { .. }, Topology::Linear) | (Ends::Circular, Topology::Circular)
            )
        {
            return Err(Error::Binding);
        }
        let stereo = d
            .links
            .iter()
            .filter(|l| {
                matches!(
                    l.linkage,
                    Linkage::Phosphorothioate {
                        stereo: Stereo::Unspecified
                    }
                )
            })
            .map(|l| l.after)
            .collect::<Vec<_>>();
        if stereo != self.unspecified_stereocenters {
            return Err(Error::Binding);
        }
        let expected_unreviewed = self
            .retained_evidence
            .iter()
            .filter(|e| match e {
                RetainedEvidence::Feature { index } => !reviews.contains(index),
                RetainedEvidence::Metadata { .. } => true,
            })
            .cloned()
            .collect::<Vec<_>>();
        if expected_unreviewed != self.unreviewed_evidence {
            return Err(Error::Binding);
        }
        if self.formula_complete
            && (d.coverage != Coverage::ClosedModel
                || !self.unreviewed_evidence.is_empty()
                || d.sites.iter().any(|s| {
                    matches!(s.nucleoside, Nucleoside::Unknown { .. })
                        || matches!(s.sugar, Sugar::Unknown { .. })
                })
                || d.links
                    .iter()
                    .any(|l| matches!(l.linkage, Linkage::Unknown { .. }))
                || matches!(
                    d.ends,
                    Ends::Linear {
                        five: FivePrime::Unknown { .. },
                        ..
                    } | Ends::Linear {
                        three: ThreePrime::Unknown { .. },
                        ..
                    }
                ))
        {
            return Err(Error::Binding);
        }
        if self
            .formula
            .is_some_and(|f| f.c == 0 || f.array().into_iter().any(|n| n > 100_000_000))
            || self.issues.len() > 8
            || self.issues.iter().any(|i| {
                i.count == 0
                    || i.positions.len() > 100
                    || i.positions.len() > i.count
                    || i.positions.iter().any(|&p| p >= self.source_length)
            })
        {
            return Err(Error::Binding);
        }
        Ok(())
    }
    /// Reconstruct every quantity from the independently supplied immutable source
    /// and declaration. A report/hash alone cannot establish chemical state.
    pub fn validate_for(&self, source: &Molecule, d: &Declaration) -> Result<()> {
        if self.engine != ENGINE
            || self.declaration.sites.len() > MAX_SITES
            || self.declaration.links.len() > MAX_LINKS
            || analyze(source, d)? != *self
        {
            return Err(Error::Binding);
        }
        Ok(())
    }
}
