//! Native anchored overlap-extension primer proposals from original DNA and
//! separately declared literal edits. No intended-product sequence is an input.
use crate::pcr_source::Template;
use crate::{
    Alphabet, Error, Location, Molecule, Sequence, Topology,
    assembly_edit_intent::Edit,
    assembly_recipe, overlap_extension,
    overlap_extension_pcr::{self, PrimerPair},
    overlap_primer_design::Fragment,
    primers::{self, DesignOptions, DesignReport},
    verification::{self, DesignPlan, PlanInput, Source, Status, VerificationReport},
};
use serde::{Deserialize, Serialize};
use std::io::Write;
pub const ENGINE_VERSION: &str = "overlap-extension-primer-design-1";
const MAX_SOURCE_BASES: usize = 2 << 20;
const MAX_SOURCE_BYTES: usize = 4 << 20;
const MAX_REPORT_BYTES: usize = 32 << 20;
const MAX_EDIT_BASES: usize = 100_000;
pub(crate) const MAX_WORK: usize = 100_000_000;
macro_rules! dto {($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Options {
    pub primer_options:DesignOptions,
    pub homology_flank:usize,
    pub max_oligo_length:usize,
    pub outer_forward_tail:String,
    pub outer_reverse_tail:String,
    pub min_paired_bases:usize,
    pub max_comparisons:usize,
    pub max_candidates:usize
}}
impl Default for Options {
    fn default() -> Self {
        Self {
            primer_options: DesignOptions::default(),
            homology_flank: 12,
            max_oligo_length: primers::MAX_OLIGO,
            outer_forward_tail: String::new(),
            outer_reverse_tail: String::new(),
            min_paired_bases: 12,
            max_comparisons: 20_000_000,
            max_candidates: 100,
        }
    }
}
dto! {pub struct Request {pub fragment:Fragment,pub edits:Vec<Edit>,pub options:Options}}
dto! {pub struct PairSearch {pub call:String,pub candidates:usize,pub examined:usize,pub rejected:usize,pub selected_index:usize,pub capped:bool,pub unique_geometry:bool,pub notes:Vec<String>}}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Report {
    pub context_revision: String,
    pub engine_version: String,
    pub source_revision: String,
    pub request: Request,
    pub preparation: overlap_extension_pcr::Request,
    /// Independently edited reference in ORIGINAL source orientation, with only
    /// the explicitly declared outer tails. Apply fragment.reverse at assembly.
    pub expected: Molecule,
    pub plan: DesignPlan,
    pub verification: VerificationReport,
    /// Three native calls on the same original molecule, not cloned source DNA.
    pub designs: Vec<DesignReport>,
    pub search: Vec<PairSearch>,
    pub warnings: Vec<String>,
}
/// Genuine native primer-planning/compilation phase, without final molecular
/// verification. Embedding workflows supply their independently frozen target
/// and verify their composed plan once. This is not physical source authority.
#[derive(Clone, Debug)]
pub struct Compilation {
    pub context_revision: String,
    pub engine_version: String,
    pub source_revision: String,
    pub request: Request,
    pub preparation: overlap_extension_pcr::Request,
    pub expected: Molecule,
    pub designs: Vec<DesignReport>,
    pub search: Vec<PairSearch>,
    pub warnings: Vec<String>,
}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("overlap-extension primer planning: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn token(s: &str) -> Result<(), Error> {
    if s.trim().is_empty() || s.len() > 128 || s.chars().any(char::is_control) {
        Err(bad("declare bounded contribution and input identities"))
    } else {
        Ok(())
    }
}
fn concrete(s: &str) -> bool {
    s.bytes().all(|b| b"ACGT".contains(&b))
}
struct Count {
    bytes: usize,
    maximum: usize,
}
impl Write for Count {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.bytes = self
            .bytes
            .checked_add(b.len())
            .filter(|n| *n <= self.maximum)
            .ok_or_else(|| std::io::Error::other("bounded overlap-extension primer report"))?;
        Ok(b.len())
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}
fn bytes(value: &impl Serialize, maximum: usize) -> Result<(), Error> {
    serde_json::to_writer(Count { bytes: 0, maximum }, value).map_err(|_| Error::TooLarge)
}
fn options(o: &Options) -> Result<(), Error> {
    let p = &o.primer_options;
    if p.min_length < 8
        || p.max_length > 80
        || p.min_length > p.max_length
        || p.max_length - p.min_length > 63
        || p.max_pairs == 0
        || p.max_pairs > 20
        || !p.forward_tail.is_empty()
        || !p.reverse_tail.is_empty()
        || !p.target_tm_c.is_finite()
        || !(-50.0..=100.0).contains(&p.target_tm_c)
        || !p.max_tm_difference_c.is_finite()
        || !(0.0..=100.0).contains(&p.max_tm_difference_c)
        || !p.min_gc_fraction.is_finite()
        || !p.max_gc_fraction.is_finite()
        || !(0.0..=1.0).contains(&p.min_gc_fraction)
        || !(p.min_gc_fraction..=1.0).contains(&p.max_gc_fraction)
        || p.max_homopolymer == 0
        || p.max_homopolymer > primers::MAX_OLIGO
        || p.max_three_prime_complementarity > primers::MAX_OLIGO
        || !p.conditions.sodium_molar.is_finite()
        || !(0.01..=1.0).contains(&p.conditions.sodium_molar)
        || !p.conditions.total_strand_molar.is_finite()
        || !(1e-12..=1e-2).contains(&p.conditions.total_strand_molar)
    {
        return Err(bad(
            "invalid native length/Tm/GC/complementarity/salt filters; base primer tails must be empty and max_pairs1..20",
        ));
    }
    if o.homology_flank == 0
        || o.homology_flank > primers::MAX_OLIGO
        || o.max_oligo_length < p.min_length
        || o.max_oligo_length > primers::MAX_OLIGO
        || o.outer_forward_tail.len() > o.max_oligo_length
        || o.outer_reverse_tail.len() > o.max_oligo_length
        || !concrete(&o.outer_forward_tail)
        || !concrete(&o.outer_reverse_tail)
        || o.min_paired_bases == 0
        || o.max_comparisons == 0
        || o.max_comparisons > overlap_extension::MAX_COMPARISONS
        || o.max_candidates == 0
        || o.max_candidates > overlap_extension::MAX_CANDIDATES
    {
        return Err(bad(
            "invalid homology, concrete outer tails, total oligo length or bounded terminal search settings",
        ));
    }
    Ok(())
}
fn conflicts(a: usize, b: usize, c: usize, d: usize) -> bool {
    if a == b {
        c <= a && a <= d
    } else if c == d {
        a <= c && c <= b
    } else {
        a < d && c < b
    }
}
struct StoredEdit<'a> {
    start: usize,
    end: usize,
    edit: &'a Edit,
}
struct Bounds<'a> {
    length: usize,
    start: usize,
    end: usize,
    insert_length: usize,
    homology: usize,
    edits: Vec<StoredEdit<'a>>,
}
fn bounds<'a>(source: Template<'_>, r: &'a Request) -> Result<Bounds<'a>, Error> {
    options(&r.options)?;
    token(&r.fragment.role)?;
    token(&r.fragment.input_id)?;
    if r.edits.is_empty() || r.edits.len() > crate::assembly_edit_intent::MAX_EDITS {
        return Err(bad("declare1..256 literal contribution edits"));
    }
    crate::model::bounds::model_preflight(&[source.reference])?;
    if source.sequence.alphabet() != Alphabet::Dna
        || source.sequence.is_empty()
        || source.sequence.len() > MAX_SOURCE_BASES
    {
        return Err(Error::TooLarge);
    }
    source.reference.validate()?;
    bytes(
        &(source.reference, source.sequence, source.strand),
        MAX_SOURCE_BYTES,
    )?;
    let f = &r.fragment;
    let n = source.sequence.len();
    let length = f
        .end
        .checked_sub(f.start)
        .filter(|n| *n > 0 && *n <= source.sequence.len())
        .ok_or_else(|| bad("select one nonempty source contribution of at most one turn"))?;
    if f.start >= n
        || f.end > f.start.checked_add(n).ok_or(Error::TooLarge)?
        || source.topology == Topology::Linear && f.end > n
    {
        return Err(bad(
            "only circular source intervals may cross the original origin",
        ));
    }
    if (0..length)
        .any(|j| !b"ACGT".contains(&source.sequence.as_str().as_bytes()[(f.start + j) % n]))
    {
        return Err(bad(
            "the selected contribution must contain concrete DNA calls",
        ));
    }
    let mut edits = Vec::with_capacity(r.edits.len());
    let mut total = 0usize;
    for (i, e) in r.edits.iter().enumerate() {
        if e.role != f.role
            || e.end.checked_sub(e.start) != Some(e.expected.len())
            || e.end > length
            || !concrete(&e.expected)
            || !concrete(&e.replacement)
            || e.expected == e.replacement
        {
            return Err(bad(
                "edits require this contribution role, bounded original half-open coordinates and distinct uppercase DNA literals",
            ));
        }
        total = total
            .checked_add(e.expected.len())
            .and_then(|n| n.checked_add(e.replacement.len()))
            .filter(|n| *n <= MAX_EDIT_BASES)
            .ok_or(Error::TooLarge)?;
        if r.edits[..i]
            .iter()
            .any(|other| conflicts(e.start, e.end, other.start, other.end))
        {
            return Err(bad(
                "edits overlap or share an ambiguous insertion boundary",
            ));
        }
        for (j, expected) in e.expected.bytes().enumerate() {
            let offset = if f.reverse {
                length - 1 - (e.start + j)
            } else {
                e.start + j
            };
            let raw = source.sequence.as_str().as_bytes()[(f.start + offset) % n];
            let actual = if f.reverse {
                crate::model::complement(raw, Alphabet::Dna)
            } else {
                raw
            };
            if actual != expected {
                return Err(bad(
                    "an expected edit literal differs from the unchanged source contribution",
                ));
            }
        }
        let (start, end) = if f.reverse {
            (length - e.end, length - e.start)
        } else {
            (e.start, e.end)
        };
        edits.push(StoredEdit {
            start,
            end,
            edit: e,
        });
    }
    edits.sort_by_key(|e| (e.start, e.end));
    let start = edits.first().ok_or(Error::TooLarge)?.start;
    let end = edits.iter().map(|e| e.end).max().ok_or(Error::TooLarge)?;
    let removed = edits.iter().try_fold(0usize, |n, e| {
        n.checked_add(e.end - e.start).ok_or(Error::TooLarge)
    })?;
    let added = edits.iter().try_fold(0usize, |n, e| {
        n.checked_add(e.edit.replacement.len())
            .ok_or(Error::TooLarge)
    })?;
    let insert_length = (end - start)
        .checked_sub(removed)
        .and_then(|n| n.checked_add(added))
        .ok_or(Error::TooLarge)?;
    let h = r.options.homology_flank;
    let minimum = r.options.primer_options.min_length;
    let right = length - end;
    if start < minimum.checked_mul(2).ok_or(Error::TooLarge)?
        || right < 2 * minimum
        || h >= start
        || h >= right
    {
        return Err(bad(
            "both unchanged flanks must fit nonoverlapping anchored primer domains and retain sequence beyond the homology; boundary edits require another preparation method",
        ));
    }
    let inner_tail = h.checked_add(insert_length).ok_or(Error::TooLarge)?;
    if inner_tail
        .checked_add(minimum)
        .is_none_or(|n| n > r.options.max_oligo_length)
        || r.options
            .outer_forward_tail
            .len()
            .checked_add(minimum)
            .is_none_or(|n| n > r.options.max_oligo_length)
        || r.options
            .outer_reverse_tail
            .len()
            .checked_add(minimum)
            .is_none_or(|n| n > r.options.max_oligo_length)
    {
        return Err(bad(
            "the central edit bridge and declared tails do not fit the oligo limit with the required annealing domains",
        ));
    }
    let homology = inner_tail.checked_add(h).ok_or(Error::TooLarge)?;
    if homology < r.options.min_paired_bases || homology > r.options.max_comparisons {
        return Err(bad(
            "selected terminal homology cannot meet the declared pairing/comparison limits",
        ));
    }
    let bounds = Bounds {
        length,
        start,
        end,
        insert_length,
        homology,
        edits,
    };
    let derived = derived_options(source, r, &bounds, "", "", false)?;
    let reserve = reserved_work(source, r, &bounds, &derived)?;
    let minimum_search = source
        .sequence
        .len()
        .checked_mul(r.options.primer_options.min_length)
        .and_then(|n| n.checked_mul(24))
        .ok_or(Error::TooLarge)?;
    reserve
        .checked_add(minimum_search)
        .filter(|n| *n <= MAX_WORK)
        .ok_or(Error::TooLarge)?;
    Ok(bounds)
}
fn reserved_work(
    source: Template<'_>,
    r: &Request,
    b: &Bounds<'_>,
    derived: &[DesignOptions; 3],
) -> Result<usize, Error> {
    let inner = r.options.homology_flank + b.insert_length;
    let mut work = 0usize;
    for o in derived {
        let count = o.max_length - o.min_length + 1;
        let width = o.max_length + o.forward_tail.len().max(o.reverse_tail.len()).max(inner);
        let pairs = count
            .checked_mul(count)
            .and_then(|n| n.checked_mul(width))
            .and_then(|n| n.checked_mul(width))
            .and_then(|n| n.checked_mul(2))
            .ok_or(Error::TooLarge)?;
        if pairs > 20_000_000 {
            return Err(bad(
                "native primer-pair work exceeds its comparison bound; narrow the length range",
            ));
        }
        let screens = count
            .checked_mul(width)
            .and_then(|n| n.checked_mul(width))
            .and_then(|n| n.checked_mul(12))
            .ok_or(Error::TooLarge)?;
        work = work
            .checked_add(pairs)
            .and_then(|n| n.checked_add(screens))
            .filter(|n| *n <= MAX_WORK)
            .ok_or(Error::TooLarge)?;
    }
    // Reserve all final verifier PCR scans and two terminal scans (planning and
    // final reconstruction). Ranked-pair scans must fit the remaining budget.
    let final_scan = source
        .comparison_bases()?
        .checked_add(b.insert_length)
        .and_then(|n| {
            n.checked_add(r.options.outer_forward_tail.len() + r.options.outer_reverse_tail.len())
        })
        .and_then(|n| n.checked_mul(r.options.primer_options.max_length))
        .and_then(|n| n.checked_mul(32))
        .ok_or(Error::TooLarge)?;
    let terminal = b
        .length
        .checked_mul(6)
        .and_then(|n| n.checked_add(b.homology))
        .ok_or(Error::TooLarge)?
        .min(r.options.max_comparisons);
    work.checked_add(final_scan)
        .and_then(|n| n.checked_add(terminal.checked_mul(2)?))
        .filter(|n| *n <= MAX_WORK)
        .ok_or(Error::TooLarge)
}
/// Borrowed source, literal-edit, geometry and aggregate work guards before any
/// full source is cloned or native primer/PCR search is started.
pub fn preflight(source: &Molecule, request: &Request) -> Result<(), Error> {
    bounds(Template::reference(source)?, request).map(|_| ())
}
pub fn design_revision(source: &Molecule, request: &Request) -> Result<String, Error> {
    preflight(source, request)?;
    assembly_recipe::hash(&(
        ENGINE_VERSION,
        Template::reference(source)?.context_revision()?,
        request,
    ))
}
fn segment(source: Template<'_>, f: &Fragment, start: usize, end: usize, out: &mut String) {
    let n = source.sequence.len();
    for i in start..end {
        out.push(source.sequence.as_str().as_bytes()[(f.start + i) % n] as char);
    }
}
fn replacement(out: &mut String, e: &Edit, reverse: bool) {
    if reverse {
        for b in e.replacement.bytes().rev() {
            out.push(crate::model::complement(b, Alphabet::Dna) as char);
        }
    } else {
        out.push_str(&e.replacement);
    }
}
fn edited(source: Template<'_>, r: &Request, b: &Bounds<'_>, start: usize, end: usize) -> String {
    let mut out = String::new();
    let mut cursor = start;
    for e in &b.edits {
        segment(source, &r.fragment, cursor, e.start, &mut out);
        replacement(&mut out, e.edit, r.fragment.reverse);
        cursor = e.end;
    }
    segment(source, &r.fragment, cursor, end, &mut out);
    out
}
fn expected(source: Template<'_>, r: &Request, b: &Bounds<'_>) -> Result<Molecule, Error> {
    let mut dna = r.options.outer_forward_tail.clone();
    dna.push_str(&edited(source, r, b, 0, b.length));
    // An absent outer tail contributes no bases; it is not a standalone DNA
    // sequence. Sequence::dna deliberately rejects an empty sequence.
    if !r.options.outer_reverse_tail.is_empty() {
        dna.push_str(
            Sequence::dna(&r.options.outer_reverse_tail)?
                .reverse_complement()?
                .as_str(),
        );
    }
    Molecule::dna(
        format!("{} independently edited preparation", r.fragment.role),
        &dna,
        Topology::Linear,
    )
}
fn derived_options(
    source: Template<'_>,
    r: &Request,
    b: &Bounds<'_>,
    left_tail: &str,
    right_tail: &str,
    real: bool,
) -> Result<[DesignOptions; 3], Error> {
    let n = source.sequence.len();
    let o = &r.options;
    let inner = o.homology_flank + b.insert_length;
    let mut left = o.primer_options.clone();
    let mut right = left.clone();
    let mut outer = left.clone();
    left.target_start = r.fragment.start;
    left.target_end = r.fragment.start + b.start;
    right.target_start = (r.fragment.start + b.end) % n;
    right.target_end = right.target_start + b.length - b.end;
    outer.target_start = r.fragment.start;
    outer.target_end = r.fragment.end;
    left.max_length = left
        .max_length
        .min(b.start - left.min_length)
        .min(o.max_oligo_length - inner);
    right.max_length = right
        .max_length
        .min(b.length - b.end - right.min_length)
        .min(o.max_oligo_length - inner);
    outer.max_length = outer
        .max_length
        .min(b.start)
        .min(b.length - b.end)
        .min(o.max_oligo_length - o.outer_forward_tail.len().max(o.outer_reverse_tail.len()));
    left.reverse_tail = if real {
        left_tail.into()
    } else {
        "A".repeat(inner)
    };
    right.forward_tail = if real {
        right_tail.into()
    } else {
        "A".repeat(inner)
    };
    outer.forward_tail = o.outer_forward_tail.clone();
    outer.reverse_tail = o.outer_reverse_tail.clone();
    for o in [&mut left, &mut right, &mut outer] {
        o.max_pairs = r.options.primer_options.max_pairs;
        if o.max_length < o.min_length {
            return Err(bad(
                "unchanged endpoints cannot fit the requested primer length filters",
            ));
        }
    }
    Ok([left, right, outer])
}
pub(crate) fn propose(
    source: Template<'_>,
    options: &DesignOptions,
    label: &str,
    role: &str,
) -> Result<DesignReport, Error> {
    let mut report =
        primers::design_primers_view(source.view(), options, &source.reference.revision())?;
    if report.pairs.is_empty() {
        return Err(bad(&format!(
            "no {label} pair passes the declared native primer filters; adjust the filters or preparation method explicitly"
        )));
    }
    for (index, pair) in report.pairs.iter_mut().enumerate() {
        for (candidate, direction) in [
            (&mut pair.forward, "forward"),
            (&mut pair.reverse, "reverse"),
        ] {
            candidate.primer.name = format!("{role}-{label}-{}-{direction}", index + 1);
            candidate.analysis.name = candidate.primer.name.clone();
        }
    }
    Ok(report)
}
pub(crate) struct Budget {
    pub(crate) used: usize,
    pub(crate) limit: usize,
    pub(crate) materialized: usize,
}
impl Budget {
    pub(crate) fn reserve(&mut self, n: usize) -> bool {
        let Some(total) = self.used.checked_add(n).filter(|n| *n <= self.limit) else {
            return false;
        };
        self.used = total;
        true
    }
    pub(crate) fn reserve_material(&mut self, n: usize) -> bool {
        let Some(total) = self
            .materialized
            .checked_add(n)
            .filter(|n| *n <= 64 * 1024 * 1024)
        else {
            return false;
        };
        self.materialized = total;
        true
    }
}
pub(crate) fn select(
    source: Template<'_>,
    d: &DesignReport,
    call: &str,
    start: usize,
    length: usize,
    budget: &mut Budget,
) -> (PairSearch, Option<crate::analysis::PcrGeometryProduct>) {
    let mut search = PairSearch {
        call: call.into(),
        candidates: d.pairs.len(),
        examined: 0,
        rejected: 0,
        selected_index: 0,
        capped: false,
        unique_geometry: false,
        notes: vec![],
    };
    for (index, pair) in d.pairs.iter().enumerate() {
        let cost = source
            .comparison_bases()
            .ok()
            .and_then(|bases| {
                bases.checked_mul(
                    pair.forward.primer.annealing.len() + pair.reverse.primer.annealing.len(),
                )
            })
            .and_then(|n| n.checked_mul(4));
        if cost.is_none_or(|cost| !budget.reserve(cost)) || budget.materialized > 32 * 1024 * 1024 {
            search.capped = true;
            search.notes.push("Aggregate source-binding comparison budget stopped ranked-pair search; the diagnostic fallback is not a uniquely established PCR.".into());
            break;
        }
        search.examined += 1;
        match source.geometries_with_usage(&pair.forward.primer, &pair.reverse.primer) {
            Ok((mut products, emitted)) => {
                budget.materialized = budget.materialized.saturating_add(emitted);
                let unique = products.len() == 1
                    && products.first().is_some_and(|g| {
                        g.product.forward_start == start % source.sequence.len()
                            && g.product.reverse_start
                                == (start + length - pair.reverse.primer.annealing.len())
                                    % source.sequence.len()
                            && g.template_length == length
                            && g.forward_tail_length == pair.forward.primer.tail.len()
                            && g.reverse_tail_length == pair.reverse.primer.tail.len()
                            && g.product.forward_primer == pair.forward.primer.name
                            && g.product.reverse_primer == pair.reverse.primer.name
                    });
                if unique {
                    search.selected_index = index;
                    search.unique_geometry = true;
                    search.notes.push(format!("Rank{} is the first searched pair with exactly one native PCR geometry matching the declared source endpoints and oligo roles. No intended DNA or target comparison selected it.",index+1));
                    return (search, products.pop());
                }
                search.rejected += 1;
                search.notes.push(format!("Rank{} retained{} PCR geometries and did not establish the unique declared source interval/oligo roles.",index+1,products.len()));
            }
            Err(error) => {
                // The native engine itself caps one call at32million bases. On
                // interruption charge that whole cap, rather than assuming it
                // allocated no competing products.
                budget.materialized = budget.materialized.saturating_add(64 * 1024 * 1024);
                search.rejected += 1;
                search.capped = true;
                search.notes.push(format!(
                    "Rank{} PCR geometry search did not complete: {}",
                    index + 1,
                    error
                ));
            }
        }
    }
    search.notes.push("No searched pair established the unique declared PCR geometry. Highest-ranked pair0 is retained as a diagnostic proposal; no product or physical authority is inferred.".into());
    (search, None)
}
fn unsearched_outer(d: &DesignReport) -> PairSearch {
    PairSearch{call:"outer".into(),candidates:d.pairs.len(),examined:0,rejected:0,selected_index:0,capped:true,unique_geometry:false,notes:vec!["Outer geometry was not searched because no unique pair of source-bound branch PCRs supplied an actual extension substrate. Highest-ranked outer pair remains diagnostic.".into()]}
}
fn preparation(
    designs: &[DesignReport],
    search: &[PairSearch],
    r: &Request,
    b: &Bounds<'_>,
    template_strand: crate::chemistry::Strand,
) -> Result<overlap_extension_pcr::Request, Error> {
    let pair = |i: usize| {
        let p = designs
            .get(i)
            .and_then(|d| d.pairs.get(search.get(i)?.selected_index))
            .ok_or_else(|| bad("missing native primer pair"))?;
        Ok::<_, Error>(PrimerPair {
            forward: p.forward.primer.clone(),
            reverse: p.reverse.primer.clone(),
            template_strand: if i == 2 {
                crate::chemistry::Strand::Top
            } else {
                template_strand
            },
        })
    };
    Ok(overlap_extension_pcr::Request {
        left: pair(0)?,
        right: pair(1)?,
        outer: pair(2)?,
        extension: overlap_extension::Options {
            overlap: b.homology,
            min_paired_bases: r.options.min_paired_bases,
            max_comparisons: r.options.max_comparisons,
            max_candidates: r.options.max_candidates,
        },
    })
}
fn plan(
    source: Template<'_>,
    r: &Request,
    expected: &Molecule,
    preparation: &overlap_extension_pcr::Request,
) -> Result<DesignPlan, Error> {
    let mut steps = vec![];
    preparation.append(Source::Input(r.fragment.input_id.clone()), &mut steps)?;
    Ok(DesignPlan {
        name: format!("{} overlap-extension PCR", r.fragment.role),
        target_revision: expected.revision(),
        inputs: vec![PlanInput {
            id: r.fragment.input_id.clone(),
            revision: source.reference.revision(),
            molecule: source.reference.clone(),
        }],
        steps,
    })
}
pub(crate) fn work_preflight(
    source: &Molecule,
    request: &Request,
) -> Result<(usize, usize), Error> {
    source_work_preflight(Template::reference(source)?, request)
}
pub(crate) fn source_work_preflight(
    source: Template<'_>,
    request: &Request,
) -> Result<(usize, usize), Error> {
    let b = bounds(source, request)?;
    let derived = derived_options(source, request, &b, "", "", false)?;
    let reserve = reserved_work(source, request, &b, &derived)?;
    let material = source
        .sequence
        .len()
        .checked_mul(6)
        .and_then(|n| n.checked_add(b.insert_length.checked_mul(8)?))
        .and_then(|n| {
            n.checked_add(
                request.options.outer_forward_tail.len() + request.options.outer_reverse_tail.len(),
            )
        })
        .ok_or(Error::TooLarge)?;
    Ok((reserve, material))
}
pub fn design(source: &Molecule, request: &Request) -> Result<Report, Error> {
    let mut budget = Budget {
        used: 0,
        limit: MAX_WORK,
        materialized: 0,
    };
    design_with_budget(source, request, &mut budget)
}
/// Shared whole-design ledger; only native source-bound proposal work uses it.
/// Plan all native primer pairs and exact operation declarations without
/// executing final verification. Source-bound ranked PCR searches still run.
pub fn compile(source: &Molecule, request: &Request) -> Result<Compilation, Error> {
    let mut budget = Budget {
        used: 0,
        limit: MAX_WORK,
        materialized: 0,
    };
    compile_source_with_budget(Template::reference(source)?, request, &mut budget)
}
pub(crate) fn design_with_budget(
    source: &Molecule,
    request: &Request,
    budget: &mut Budget,
) -> Result<Report, Error> {
    let compiled = compile_source_with_budget(Template::reference(source)?, request, budget)?;
    let actual_plan = plan(
        Template::reference(source)?,
        request,
        &compiled.expected,
        &compiled.preparation,
    )?;
    let verification = verification::verify(&actual_plan, &compiled.expected);
    let mut warnings = compiled.warnings;
    if verification.status != Status::ConsistentUnderAssumptions {
        warnings.push("The selected oligos did not establish a complete target-consistent physical preparation. Keep the verifier findings; no intended sequence was substituted to obtain success.".into());
    }
    let report = Report {
        context_revision: compiled.context_revision,
        engine_version: compiled.engine_version,
        source_revision: compiled.source_revision,
        request: compiled.request,
        preparation: compiled.preparation,
        expected: compiled.expected,
        plan: actual_plan,
        verification,
        designs: compiled.designs,
        search: compiled.search,
        warnings,
    };
    report.validate()?;
    Ok(report)
}
pub(crate) fn compile_source_with_budget(
    source: Template<'_>,
    request: &Request,
    budget: &mut Budget,
) -> Result<Compilation, Error> {
    let b = bounds(source, request)?;
    let source_revision = source.reference.revision();
    let context_revision =
        assembly_recipe::hash(&(ENGINE_VERSION, source.context_revision()?, request))?;
    let insert = edited(source, request, &b, b.start, b.end);
    let mut right_prefix = String::with_capacity(request.options.homology_flank);
    segment(
        source,
        &request.fragment,
        b.end,
        b.end + request.options.homology_flank,
        &mut right_prefix,
    );
    let left_tail = Sequence::dna(&format!("{insert}{right_prefix}"))?
        .reverse_complement()?
        .as_str()
        .to_owned();
    let mut right_tail = String::with_capacity(request.options.homology_flank + insert.len());
    segment(
        source,
        &request.fragment,
        b.start - request.options.homology_flank,
        b.start,
        &mut right_tail,
    );
    right_tail.push_str(&insert);
    let derived = derived_options(source, request, &b, &left_tail, &right_tail, true)?;
    let reserve = reserved_work(source, request, &b, &derived)?;
    let material = source
        .sequence
        .len()
        .checked_mul(6)
        .and_then(|n| n.checked_add(b.insert_length.checked_mul(8)?))
        .and_then(|n| {
            n.checked_add(
                request.options.outer_forward_tail.len() + request.options.outer_reverse_tail.len(),
            )
        })
        .ok_or(Error::TooLarge)?;
    if !budget.reserve(reserve) || !budget.reserve_material(material) {
        return Err(Error::TooLarge);
    }
    let mut designs = Vec::with_capacity(3);
    for (o, label) in derived.iter().zip(["left", "right", "outer"]) {
        designs.push(propose(source, o, label, &request.fragment.role)?);
    }
    let (left, left_product) = select(
        source,
        &designs[0],
        "left",
        derived[0].target_start,
        derived[0].target_end - derived[0].target_start,
        budget,
    );
    let (right, right_product) = select(
        source,
        &designs[1],
        "right",
        derived[1].target_start,
        derived[1].target_end - derived[1].target_start,
        budget,
    );
    let outer = if let (Some(left), Some(right)) = (left_product, right_product) {
        let a = crate::chemistry::Duplex::exact_pcr(&left.product.sequence)?;
        let z = crate::chemistry::Duplex::exact_pcr(&right.product.sequence)?;
        let o = overlap_extension::Options {
            overlap: b.homology,
            min_paired_bases: request.options.min_paired_bases,
            max_comparisons: request.options.max_comparisons,
            max_candidates: request.options.max_candidates,
        };
        match overlap_extension::geometry(&a, &z, &o) {
            Ok(geometry) => {
                let actual = geometry
                    .fill
                    .duplex
                    .molecule("Source-bound selected reciprocal extension")?;
                let tracked = crate::modified_duplex::TrackedDuplex::from_unknown_duplex(
                    &geometry.fill.duplex,
                )?;
                select(
                    Template::fresh(&actual, &tracked, crate::chemistry::Strand::Top)?,
                    &designs[2],
                    "outer",
                    0,
                    actual.sequence.len(),
                    budget,
                )
                .0
            }
            Err(error) => {
                let mut report = unsearched_outer(&designs[2]);
                report.notes.push(format!(
                    "Selected branch extension was unavailable: {error}"
                ));
                report
            }
        }
    } else {
        unsearched_outer(&designs[2])
    };
    let search = vec![left, right, outer];
    // Expected DNA is a separate source/literal construction. It is never read by
    // ranked pair selection, which uses only actual source endpoints/strand roles.
    let expected = expected(source, request, &b)?;
    let preparation = preparation(&designs, &search, request, &b, source.strand)?;
    preparation.validate_shape()?;
    let mut warnings=vec![
        "Native anchored pairs are ranked under the requested filters (at most20). Bounded search chooses the first unique declared source-endpoint/oligo-role PCR geometry in each ranked list; the outer list is checked on actual extended branch material. No intended target DNA chooses a pair or product. Searches are not globally optimal or exhaustive beyond the retained native lists.".into(),
        "Derived endpoints, physical maximum lengths, tails and max_pairs are recorded per native call. Requested Tm/GC/salt/structure filters remain unchanged. Work caps or ambiguous/competing PCRs retain an explicit diagnostic fallback rather than inventing material.".into(),
        "Literal edits use original construct-oriented contribution coordinates; reverse edits and replacement DNA are converted back to stored source orientation. The prepared product remains in stored source orientation; apply fragment.reverse only at subsequent assembly.".into(),
        "Expected DNA is independently built from the unchanged source contribution, explicit expected/replacement literals and declared outer tails. It never supplies PCR template material. Reference-geometry uniqueness does not establish physical source chemistry, terminal-mixture uniqueness or experimental success.".into(),
        "Outer forward/reverse tails refer to oligos in stored source orientation. Their endpoints remain on unchanged original flanks. Readable saved-material metadata cannot restore source chemistry; exact closed-origin reconstruction remains necessary where the verifier requires it.".into(),
    ];
    warnings.push(format!("Ranked source searches charged{} comparison-work units (including reserved native design/final-verification work) within{} and{} materialized-base units within64MiB. Native interrupted materialization calls are conservatively charged at their full32MiB cap.",budget.used,MAX_WORK,budget.materialized));
    let report = Compilation {
        context_revision,
        engine_version: ENGINE_VERSION.into(),
        source_revision,
        request: request.clone(),
        preparation,
        expected,
        designs,
        search,
        warnings,
    };
    validate_source_preparation(
        source,
        request,
        &report.context_revision,
        &report.preparation,
        &report.designs,
        &report.search,
    )?;
    bytes(
        &(
            &report.context_revision,
            &report.request,
            &report.preparation,
            &report.expected,
            &report.designs,
            &report.search,
            &report.warnings,
        ),
        MAX_REPORT_BYTES,
    )?;
    Ok(report)
}
fn location_guard(location: &Location, n: usize) -> Result<(), Error> {
    let mut pending = vec![(location, 0usize)];
    let mut count = 0usize;
    while let Some((l, depth)) = pending.pop() {
        count += 1;
        if count > 8 || depth > 4 {
            return Err(Error::TooLarge);
        }
        match l {
            Location::Range {
                start,
                end,
                before: false,
                after: false,
            } if start < end && *end <= n => {}
            Location::Complement(inner) => pending.push((inner, depth + 1)),
            Location::Join(parts) if parts.len() <= 2 => {
                for p in parts {
                    pending.push((p, depth + 1));
                }
            }
            _ => return Err(bad("invalid native primer location geometry")),
        }
    }
    Ok(())
}
/// Bounded readonly correspondence of actual native calls, selected primers
/// and literal edits; no primer ranking/PCR/verification is repeated.
pub fn validate_preparation(
    source: &Molecule,
    request: &Request,
    context: &str,
    prepared: &overlap_extension_pcr::Request,
    designs: &[DesignReport],
    searches: &[PairSearch],
) -> Result<Molecule, Error> {
    validate_source_preparation(
        Template::reference(source)?,
        request,
        context,
        prepared,
        designs,
        searches,
    )
}
pub(crate) fn validate_source_preparation(
    source: Template<'_>,
    request: &Request,
    context: &str,
    prepared: &overlap_extension_pcr::Request,
    designs: &[DesignReport],
    searches: &[PairSearch],
) -> Result<Molecule, Error> {
    if designs.len() != 3 || searches.len() != 3 {
        return Err(bad("expected three native primer calls"));
    }
    for (d, search) in designs.iter().zip(searches) {
        if d.pairs.is_empty()
            || d.pairs.len() > 20
            || search.candidates != d.pairs.len()
            || search.examined > search.candidates
            || search.rejected > search.examined
            || search.selected_index >= d.pairs.len()
            || search.notes.len() > 22
            || search.notes.iter().any(|s| s.len() > 4096)
            || search.unique_geometry
                && (search.examined != search.selected_index + 1
                    || search.rejected + 1 != search.examined)
            || !search.unique_geometry
                && (search.selected_index != 0 || search.rejected != search.examined)
        {
            return Err(bad("invalid bounded ranked-pair search evidence"));
        }
    }
    bytes(&(request, prepared, designs, searches), MAX_REPORT_BYTES)?;
    let b = bounds(source, request)?;
    let source_revision = source.reference.revision();
    if assembly_recipe::hash(&(ENGINE_VERSION, source.context_revision()?, request))? != context {
        return Err(bad("returned source/request binding changed"));
    }
    prepared.validate_shape()?;
    let insert = edited(source, &request, &b, b.start, b.end);
    let mut right_prefix = String::new();
    segment(
        source,
        &request.fragment,
        b.end,
        b.end + request.options.homology_flank,
        &mut right_prefix,
    );
    let left_tail = Sequence::dna(&format!("{insert}{right_prefix}"))?
        .reverse_complement()?
        .as_str()
        .to_owned();
    let mut right_tail = String::new();
    segment(
        source,
        &request.fragment,
        b.start - request.options.homology_flank,
        b.start,
        &mut right_tail,
    );
    right_tail.push_str(&insert);
    let derived = derived_options(source, &request, &b, &left_tail, &right_tail, true)?;
    for (call_index, (d, o)) in designs.iter().zip(&derived).enumerate() {
        if d.template_strand != source.strand
            || d.source_revision != source_revision
            || d.options != *o
            || d.pairs.len() > o.max_pairs
            || searches[call_index].call != ["left", "right", "outer"][call_index]
            || d.pairs.windows(2).any(|pair| pair[0].score > pair[1].score)
            || d.evaluated_pairs > 4096
            || d.rejected_pairs > d.evaluated_pairs
            || d.warnings.len() > 32
            || d.warnings.iter().any(|s| s.len() > 8192)
        {
            return Err(bad("returned native primer-call settings changed"));
        }
        for pair in &d.pairs {
            if !pair.score.is_finite()
                || !pair.tm_difference_c.is_finite()
                || pair.forward.primer.annealing.len() + pair.reverse.primer.annealing.len()
                    > o.target_end - o.target_start
                || pair.product_length
                    != (o.target_end - o.target_start) + o.forward_tail.len() + o.reverse_tail.len()
            {
                return Err(bad("returned primer pair shape changed"));
            }
            for (candidate, reverse, tail) in [
                (&pair.forward, false, &o.forward_tail),
                (&pair.reverse, true, &o.reverse_tail),
            ] {
                let p = &candidate.primer;
                overlap_extension_pcr::primer(p)?;
                if p.annealing.len() < o.min_length
                    || p.annealing.len() > o.max_length
                    || p.tail != *tail
                    || p.ordered_sequence().len() > request.options.max_oligo_length
                    || candidate.analysis.name != p.name
                    || candidate.analysis.annealing != p.annealing.as_str()
                    || candidate.analysis.tail != p.tail
                    || candidate.analysis.ordered_sequence != p.ordered_sequence()
                {
                    return Err(bad("returned oligo differs from its native call/analysis"));
                }
                location_guard(&candidate.location, source.sequence.len())?;
                let mut actual = String::with_capacity(p.annealing.len());
                let start = if reverse {
                    o.target_end - p.annealing.len()
                } else {
                    o.target_start
                };
                for j in 0..p.annealing.len() {
                    actual.push(
                        source.sequence.as_str().as_bytes()[(start + j) % source.sequence.len()]
                            as char,
                    );
                }
                if reverse {
                    actual = Sequence::dna(&actual)?
                        .reverse_complement()?
                        .as_str()
                        .to_owned();
                }
                if p.annealing.as_str() != actual
                    || candidate.location.extract(&source.sequence)? != actual
                {
                    return Err(bad(
                        "returned primer does not bind its original anchored endpoint",
                    ));
                }
            }
        }
    }
    let selected = preparation(designs, searches, request, &b, source.strand)?;
    if assembly_recipe::hash(&selected)? != assembly_recipe::hash(prepared)? {
        return Err(bad(
            "returned preparation differs from the selected native oligos",
        ));
    }
    expected(source, request, &b)
}
impl Report {
    /// Bounded documentary correspondence only. This does not repeat primer
    /// ranking/PCR or grant material authority to an imported report.
    pub fn validate(&self) -> Result<(), Error> {
        if !hex(&self.context_revision)
            || !hex(&self.source_revision)
            || self.engine_version != ENGINE_VERSION
            || self.plan.inputs.len() != 1
            || self.plan.steps.len() != 4
            || self.designs.len() != 3
            || self.search.len() != 3
            || self.verification.products.len() > 4
            || self.verification.physical_processing_reports.len() > 1
            || self.verification.step_reports.len() > 4
            || self.verification.processed_step_reports.len() > 4
            || self.verification.reference_step_reports.len() > 4
            || self.verification.findings.len() > 1000
            || self.verification.assumptions.len() > 100
            || self
                .verification
                .findings
                .iter()
                .chain(&self.verification.assumptions)
                .any(|s| s.len() > 8192)
            || self.warnings.len() > 32
            || self.warnings.iter().any(|s| s.len() > 8192)
        {
            return Err(bad("invalid returned report bounds or version"));
        }
        bytes(self, MAX_REPORT_BYTES)?;
        let input = &self.plan.inputs[0];
        let source = &input.molecule;
        if input.id != self.request.fragment.input_id
            || input.revision != self.source_revision
            || source.revision() != self.source_revision
        {
            return Err(bad("returned source input changed"));
        }
        let expected = validate_preparation(
            source,
            &self.request,
            &self.context_revision,
            &self.preparation,
            &self.designs,
            &self.search,
        )?;
        let preparation = &self.preparation;
        if expected.revision() != self.expected.revision() {
            return Err(bad(
                "returned expected reference differs from independent literal edits",
            ));
        }
        let mut expected_steps = vec![];
        preparation.append(
            Source::Input(self.request.fragment.input_id.clone()),
            &mut expected_steps,
        )?;
        if self.plan.name != format!("{} overlap-extension PCR", self.request.fragment.role)
            || self.plan.target_revision != expected.revision()
            || assembly_recipe::hash(&expected_steps)? != assembly_recipe::hash(&self.plan.steps)?
            || self.verification.plan_revision != verification::plan_revision(&self.plan)
            || self.verification.target_revision != self.expected.revision()
            || self.verification.engine_version != verification::ENGINE_VERSION
            || self.verification.status == Status::ConsistentUnderAssumptions
                && (self.verification.products.len() != 4
                    || !self.verification.products.last().is_some_and(|p| {
                        p.sequence == self.expected.sequence && p.topology == Topology::Linear
                    }))
        {
            return Err(bad(
                "returned plan/verification/expected-reference correspondence changed",
            ));
        }
        for p in &self.verification.products {
            crate::model::bounds::model_preflight(&[p])?;
        }
        bytes(self, MAX_REPORT_BYTES)
    }
    pub fn validate_binding(&self, revision: &str) -> Result<(), Error> {
        if self.context_revision != revision {
            return Err(bad("stale source/edit/options context"));
        }
        self.validate()
    }
}
