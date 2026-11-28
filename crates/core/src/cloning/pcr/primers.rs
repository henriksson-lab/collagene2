//! Independent basic oligo calculations and anchored primer selection.
//! This is not Primer3, a folding engine or an experimental PCR-success predictor.
use crate::sequence::analysis::{self, Primer};
use collagene_model::{Alphabet, Error, Location, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};

pub const MAX_OLIGO: usize = 200;
const MAX_BINDING_WORK: usize = 200_000_000;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct ThermoConditions {
    /// Sodium concentration in mol/litre, not a magnesium-equivalent salt value.
    pub sodium_molar: f64,
    /// Sum of both strands for the equal-concentration perfect-duplex model.
    pub total_strand_molar: f64,
}
impl Default for ThermoConditions {
    fn default() -> Self {
        Self {
            sodium_molar: 0.05,
            total_strand_molar: 500e-9,
        }
    }
}
impl ThermoConditions {
    fn validate(&self) -> Result<(), Error> {
        if !self.sodium_molar.is_finite()
            || !(0.01..=1.0).contains(&self.sodium_molar)
            || !self.total_strand_molar.is_finite()
            || !(1e-12..=1e-2).contains(&self.total_strand_molar)
        {
            return Err(Error::Unsupported(
                "basic Tm requires sodium 0.01–1 M and total strand concentration 1 pM–10 mM"
                    .into(),
            ));
        }
        Ok(())
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GcBounds {
    pub minimum_fraction: f64,
    pub maximum_fraction: f64,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TmReport {
    pub model: String,
    pub conditions: ThermoConditions,
    pub tm_c: f64,
    pub delta_h_kcal_per_mol: f64,
    pub delta_s_cal_per_mol_kelvin_at_1m: f64,
    pub salt_corrected_delta_s: f64,
    pub self_complementary: bool,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComplementaryRun {
    /// Both intervals use each oligo's own 5′→3′ coordinates. Pairing is antiparallel.
    pub a_start: usize,
    pub a_end: usize,
    pub b_start: usize,
    pub b_end: usize,
    pub length: usize,
    pub certain: bool,
    pub a_three_prime: bool,
    pub b_three_prime: bool,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ComplementarityReport {
    pub max_run: usize,
    /// Longest contiguous run involving at least one physical 3′ terminus.
    pub max_three_prime_run: usize,
    pub strongest: Option<ComplementaryRun>,
    pub warnings: Vec<String>,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrimerAnalysis {
    pub name: String,
    pub annealing: String,
    pub tail: String,
    pub ordered_sequence: String,
    pub annealing_gc: GcBounds,
    pub ordered_gc: GcBounds,
    pub tm: Option<TmReport>,
    pub tm_unavailable: Option<String>,
    pub homopolymer_run: usize,
    pub hairpin: ComplementarityReport,
    pub homodimer: ComplementarityReport,
    pub warnings: Vec<String>,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct BindingOptions {
    pub max_mismatches: usize,
    pub min_three_prime_exact: usize,
    pub max_hits: usize,
    pub max_work: usize,
}
impl Default for BindingOptions {
    fn default() -> Self {
        Self {
            max_mismatches: 0,
            min_three_prime_exact: 0,
            max_hits: 1000,
            max_work: 20_000_000,
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BindingHit {
    pub start: usize,
    pub length: usize,
    pub reverse: bool,
    pub crosses_origin: bool,
    pub certain: bool,
    /// Disjoint and compatible-ambiguous bases, in the primer annealing domain's own 5′→3′ coordinates.
    pub mismatches: Vec<usize>,
    pub ambiguous_positions: Vec<usize>,
    pub three_prime_exact: usize,
    pub location: Location,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BindingReport {
    pub source_revision: String,
    pub primer_name: String,
    pub annealing: String,
    pub options: BindingOptions,
    pub hits: Vec<BindingHit>,
    pub truncated: bool,
    pub placements_examined: usize,
    pub warnings: Vec<String>,
}
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct DesignOptions {
    pub target_start: usize,
    pub target_end: usize,
    pub min_length: usize,
    pub max_length: usize,
    pub target_tm_c: f64,
    pub max_tm_difference_c: f64,
    pub min_gc_fraction: f64,
    pub max_gc_fraction: f64,
    pub max_homopolymer: usize,
    pub max_three_prime_complementarity: usize,
    pub forward_tail: String,
    pub reverse_tail: String,
    pub conditions: ThermoConditions,
    pub max_pairs: usize,
}
impl Default for DesignOptions {
    fn default() -> Self {
        Self {
            target_start: 0,
            target_end: 100,
            min_length: 18,
            max_length: 30,
            target_tm_c: 60.,
            max_tm_difference_c: 5.,
            min_gc_fraction: 0.2,
            max_gc_fraction: 0.8,
            max_homopolymer: 5,
            max_three_prime_complementarity: 3,
            forward_tail: String::new(),
            reverse_tail: String::new(),
            conditions: ThermoConditions::default(),
            max_pairs: 20,
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrimerCandidate {
    pub primer: Primer,
    pub location: Location,
    pub analysis: PrimerAnalysis,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrimerPair {
    pub forward: PrimerCandidate,
    pub reverse: PrimerCandidate,
    pub heterodimer: ComplementarityReport,
    pub tm_difference_c: f64,
    pub score: f64,
    pub product_length: usize,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DesignReport {
    pub template_strand: crate::chemistry::Strand,
    pub source_revision: String,
    pub options: DesignOptions,
    pub pairs: Vec<PrimerPair>,
    pub evaluated_pairs: usize,
    pub rejected_pairs: usize,
    pub warnings: Vec<String>,
}

fn checked_primer(primer: &Primer) -> Result<Primer, Error> {
    if primer.name.len() > 1024
        || primer.annealing.len() > MAX_OLIGO
        || primer.tail.len() > MAX_OLIGO
        || primer.annealing.len().saturating_add(primer.tail.len()) > MAX_OLIGO
    {
        return Err(Error::TooLarge);
    }
    if primer.annealing.alphabet() != Alphabet::Dna {
        return Err(Error::Unsupported(
            "oligo tools require DNA annealing domains".into(),
        ));
    }
    Primer::new(&primer.name, primer.annealing.as_str(), &primer.tail)
}
fn gc(text: &str) -> GcBounds {
    let (mut minimum, mut maximum) = (0usize, 0usize);
    for b in text.bytes() {
        let mask = analysis::base_mask(b);
        minimum += usize::from(mask & !6 == 0);
        maximum += usize::from(mask & 6 != 0);
    }
    GcBounds {
        minimum_fraction: minimum as f64 / text.len() as f64,
        maximum_fraction: maximum as f64 / text.len() as f64,
    }
}
fn concrete(text: &str) -> bool {
    text.bytes()
        .all(|b| analysis::base_mask(b).count_ones() == 1)
}
fn complement_mask(mask: u8) -> u8 {
    ((mask & 1) << 3) | ((mask & 2) << 1) | ((mask & 4) >> 1) | ((mask & 8) >> 3)
}
fn homopolymer(text: &str) -> usize {
    let (mut best, mut current, mut last) = (0, 0, 0);
    for b in text.bytes() {
        if b == last && analysis::base_mask(b).count_ones() == 1 {
            current += 1;
        } else {
            current = usize::from(analysis::base_mask(b).count_ones() == 1);
        }
        best = best.max(current);
        last = b;
    }
    best
}
fn nn(text: &str, conditions: &ThermoConditions) -> Result<TmReport, Error> {
    conditions.validate()?;
    if text.len() < 2 || !concrete(text) {
        return Err(Error::Unsupported(
            "perfect-duplex Tm requires at least two unambiguous DNA bases".into(),
        ));
    }
    // Factual values from SantaLucia1998 Table2, not copied from software.
    let (mut h, mut s) = (0., 0.);
    for end in [text.as_bytes()[0], text.as_bytes()[text.len() - 1]] {
        let (eh, es) = if end == b'A' || end == b'T' {
            (2.3, 4.1)
        } else {
            (0.1, -2.8)
        };
        h += eh;
        s += es;
    }
    for pair in text.as_bytes().windows(2) {
        let (ph, ps) = match pair {
            b"AA" | b"TT" => (-7.9, -22.2),
            b"AT" => (-7.2, -20.4),
            b"TA" => (-7.2, -21.3),
            b"CA" | b"TG" => (-8.5, -22.7),
            b"GT" | b"AC" => (-8.4, -22.4),
            b"CT" | b"AG" => (-7.8, -21.0),
            b"GA" | b"TC" => (-8.2, -22.2),
            b"CG" => (-10.6, -27.2),
            b"GC" => (-9.8, -24.4),
            b"GG" | b"CC" => (-8.0, -19.9),
            _ => unreachable!("validated concrete DNA"),
        };
        h += ph;
        s += ps;
    }
    let self_complementary = Sequence::dna(text)?.reverse_complement()?.as_str() == text;
    if self_complementary {
        s -= 1.4;
    }
    let salt_s = s + 0.368 * (text.len() - 1) as f64 * conditions.sodium_molar.ln();
    let concentration = conditions.total_strand_molar / if self_complementary { 1. } else { 4. };
    let tm = h * 1000. / (salt_s + 1.987 * concentration.ln()) - 273.15;
    if !tm.is_finite() {
        return Err(Error::Unsupported(
            "Tm is undefined under these conditions".into(),
        ));
    }
    Ok(TmReport {
        model: "SantaLucia1998 perfect DNA duplex; sodium entropy correction".into(),
        conditions: conditions.clone(),
        tm_c: tm,
        delta_h_kcal_per_mol: h,
        delta_s_cal_per_mol_kelvin_at_1m: s,
        salt_corrected_delta_s: salt_s,
        self_complementary,
    })
}
fn screen(a: &str, b: &str, hairpin: bool) -> ComplementarityReport {
    let (a, b) = (a.as_bytes(), b.as_bytes());
    let mut report=ComplementarityReport {max_run:0,max_three_prime_run:0,strongest:None,warnings:vec!["Ungapped contiguous complementarity screen only; no folding energy, loop/bulge/mismatch structure or polymerase extension probability is predicted.".into()]};
    for offset in -(b.len() as isize - 1)..a.len() as isize {
        let (mut run, mut certain) = (0usize, true);
        for i in 0..a.len() {
            let br = i as isize - offset;
            if br < 0 || br >= b.len() as isize {
                run = 0;
                certain = true;
                continue;
            }
            let j = b.len() - 1 - br as usize;
            let valid = !hairpin || (j > i && j - i > 3); // minimum three unpaired loop bases
            let am = analysis::base_mask(a[i]);
            let bm = complement_mask(analysis::base_mask(b[j]));
            if !valid || am & bm == 0 {
                run = 0;
                certain = true;
                continue;
            }
            run += 1;
            certain &= am.count_ones() == 1 && bm.count_ones() == 1;
            let a_start = i + 1 - run;
            let b_end = j + run;
            let a_three_prime = i + 1 == a.len();
            let b_three_prime = b_end == b.len();
            if a_three_prime || b_three_prime {
                report.max_three_prime_run = report.max_three_prime_run.max(run);
            }
            if run > report.max_run
                || run == report.max_run
                    && report
                        .strongest
                        .as_ref()
                        .is_some_and(|r| !r.certain && certain)
            {
                report.max_run = run;
                report.strongest = Some(ComplementaryRun {
                    a_start,
                    a_end: i + 1,
                    b_start: j,
                    b_end,
                    length: run,
                    certain,
                    a_three_prime,
                    b_three_prime,
                });
            }
        }
    }
    report
}
pub fn analyze_primer(
    primer: &Primer,
    conditions: &ThermoConditions,
) -> Result<PrimerAnalysis, Error> {
    conditions.validate()?;
    let primer = checked_primer(primer)?;
    let ordered = primer.ordered_sequence();
    let (tm, tm_unavailable) = match nn(primer.annealing.as_str(), conditions) {
        Ok(v) => (Some(v), None),
        Err(e) => (None, Some(e.to_string())),
    };
    let mut warnings=vec!["Tm describes an equal-concentration perfect DNA duplex of the annealing domain; it is not a recommended PCR annealing temperature. Magnesium, dNTPs, mismatches, modifications and additives are not modeled.".into(),"Secondary-structure screens include the ordered 5′ tail; the initial template-binding Tm excludes it.".into()];
    if !concrete(&ordered) {
        warnings.push("Ambiguous bases yield GC bounds and possible complementarity; no single mixture Tm is assigned.".into());
    }
    Ok(PrimerAnalysis {
        name: primer.name.clone(),
        annealing: primer.annealing.as_str().into(),
        tail: primer.tail.clone(),
        annealing_gc: gc(primer.annealing.as_str()),
        ordered_gc: gc(&ordered),
        tm,
        tm_unavailable,
        homopolymer_run: homopolymer(&ordered),
        hairpin: screen(&ordered, &ordered, true),
        homodimer: screen(&ordered, &ordered, false),
        ordered_sequence: ordered,
        warnings,
    })
}
pub fn heterodimer(a: &Primer, b: &Primer) -> Result<ComplementarityReport, Error> {
    let (a, b) = (checked_primer(a)?, checked_primer(b)?);
    Ok(screen(&a.ordered_sequence(), &b.ordered_sequence(), false))
}
fn mapped_location(start: usize, length: usize, n: usize, reverse: bool) -> Location {
    let location = if start + length <= n {
        Location::range(start, start + length)
    } else {
        Location::Join(vec![
            Location::range(start, n),
            Location::range(0, start + length - n),
        ])
    };
    if reverse {
        Location::Complement(Box::new(location))
    } else {
        location
    }
}
pub fn find_binding_sites(
    molecule: &Molecule,
    primer: &Primer,
    options: &BindingOptions,
) -> Result<BindingReport, Error> {
    molecule.validate()?;
    let primer = checked_primer(primer)?;
    if molecule.sequence.alphabet() != Alphabet::Dna {
        return Err(Error::Unsupported(
            "primer binding requires a DNA reference".into(),
        ));
    }
    let m = primer.annealing.len();
    let n = molecule.sequence.len();
    if options.max_mismatches > m
        || options.min_three_prime_exact > m
        || options.max_hits == 0
        || options.max_hits > 10_000
        || options.max_work == 0
        || options.max_work > MAX_BINDING_WORK
    {
        return Err(Error::Unsupported(
            "invalid binding mismatch/3′/hit/work limits".into(),
        ));
    }
    let mut report=BindingReport {source_revision:molecule.revision(),primer_name:primer.name.clone(),annealing:primer.annealing.as_str().into(),options:options.clone(),hits:vec![],truncated:false,placements_examined:0,warnings:vec!["Reference-sequence compatibility only: mismatch counts do not predict amplification, specificity outside this document, or experimental success. Tails are excluded from initial binding.".into()]};
    if m > n {
        return Ok(report);
    }
    let limit = if molecule.topology == Topology::Circular {
        n
    } else {
        n - m + 1
    };
    let reverse = primer.annealing.reverse_complement()?;
    let mut work = 0usize;
    'orientations: for (query, is_reverse) in
        [(primer.annealing.as_str(), false), (reverse.as_str(), true)]
    {
        for start in 0..limit {
            if work
                .checked_add(m)
                .is_none_or(|next| next > options.max_work)
            {
                report.truncated = true;
                break 'orientations;
            }
            work += m;
            report.placements_examined += 1;
            let (mut mismatches, mut ambiguous) = (Vec::new(), Vec::new());
            let mut exact = vec![false; m];
            for (j, q) in query.bytes().enumerate() {
                let a = analysis::base_mask(molecule.sequence.as_str().as_bytes()[(start + j) % n]);
                let b = analysis::base_mask(q);
                let p = if is_reverse { m - 1 - j } else { j };
                if a & b == 0 {
                    mismatches.push(p);
                } else if a.count_ones() != 1 || b.count_ones() != 1 {
                    ambiguous.push(p);
                } else {
                    exact[p] = true;
                }
            }
            let three_prime_exact = exact.iter().rev().take_while(|v| **v).count();
            if mismatches.len() > options.max_mismatches
                || three_prime_exact < options.min_three_prime_exact
            {
                continue;
            }
            if report.hits.len() == options.max_hits {
                report.truncated = true;
                break 'orientations;
            }
            mismatches.sort_unstable();
            ambiguous.sort_unstable();
            report.hits.push(BindingHit {
                start,
                length: m,
                reverse: is_reverse,
                crosses_origin: start + m > n,
                certain: mismatches.is_empty() && ambiguous.is_empty(),
                mismatches,
                ambiguous_positions: ambiguous,
                three_prime_exact,
                location: mapped_location(start, m, n, is_reverse),
            });
        }
    }
    report.hits.sort_by_key(|h| (h.start, h.reverse));
    if report.truncated {
        report.warnings.push("Binding enumeration was truncated by the work/hit budget; absence of additional sites is not established.".into());
    }
    Ok(report)
}
pub fn design_primers(molecule: &Molecule, options: &DesignOptions) -> Result<DesignReport, Error> {
    molecule.validate()?;
    design_primers_view(
        crate::sequence::analysis::TemplateView::reference(molecule),
        options,
        &molecule.revision(),
    )
}
/// Rank oligos against a borrowed actual source frame. The supplied immutable
/// source revision is documentary identity, never a replacement DNA template.
pub fn design_primers_view(
    molecule: crate::sequence::analysis::TemplateView<'_>,
    options: &DesignOptions,
    source_revision: &str,
) -> Result<DesignReport, Error> {
    if source_revision.len() != 64
        || !source_revision
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
        || molecule.sequence.len() > collagene_model::MAX_SEQUENCE_LEN
    {
        return Err(Error::TooLarge);
    }
    options.conditions.validate()?;
    if molecule.sequence.alphabet() != Alphabet::Dna {
        return Err(Error::Unsupported("primer selection requires DNA".into()));
    }
    let n = molecule.sequence.len();
    let target_len = options
        .target_end
        .checked_sub(options.target_start)
        .ok_or_else(|| Error::Coordinates("target end precedes start".into()))?;
    if options.target_start >= n
        || target_len == 0
        || target_len > n
        || molecule.topology == Topology::Linear && options.target_end > n
    {
        return Err(Error::Coordinates(
            "target interval must contain at most one source turn".into(),
        ));
    }
    if options.min_length < 8
        || options.max_length > 80
        || options.min_length > options.max_length
        || options.max_length - options.min_length > 63
        || options.max_pairs == 0
        || options.max_pairs > 100
        || !options.target_tm_c.is_finite()
        || !(-50.0..=100.0).contains(&options.target_tm_c)
        || !options.max_tm_difference_c.is_finite()
        || !(0.0..=100.0).contains(&options.max_tm_difference_c)
        || !options.min_gc_fraction.is_finite()
        || !options.max_gc_fraction.is_finite()
        || !(0.0..=1.0).contains(&options.min_gc_fraction)
        || !(options.min_gc_fraction..=1.0).contains(&options.max_gc_fraction)
        || options.max_homopolymer == 0
        || options.max_homopolymer > MAX_OLIGO
        || options.max_three_prime_complementarity > MAX_OLIGO
    {
        return Err(Error::Unsupported(
            "invalid primer selection length, Tm, GC, structure or output limits".into(),
        ));
    }
    // Validate even when no candidate lengths fit the target.
    checked_primer(&Primer::new("forward", "ACGTACGT", &options.forward_tail)?)?;
    checked_primer(&Primer::new("reverse", "ACGTACGT", &options.reverse_tail)?)?;
    let (mut forward, mut reverse) = (Vec::new(), Vec::new());
    for length in options.min_length..=options.max_length.min(target_len) {
        for (is_reverse, tail, list) in [
            (false, &options.forward_tail, &mut forward),
            (true, &options.reverse_tail, &mut reverse),
        ] {
            let start = if is_reverse {
                (options.target_end - length) % n
            } else {
                options.target_start
            };
            let location = mapped_location(start, length, n, is_reverse);
            let annealing = location.extract(&molecule.sequence)?;
            let primer = Primer::new(
                if is_reverse { "reverse" } else { "forward" },
                &annealing,
                tail,
            )?;
            let metrics = analyze_primer(&primer, &options.conditions)?;
            let gc = &metrics.annealing_gc;
            if metrics.tm.is_none()
                || gc.minimum_fraction < options.min_gc_fraction
                || gc.maximum_fraction > options.max_gc_fraction
                || metrics.homopolymer_run > options.max_homopolymer
                || metrics.homodimer.max_three_prime_run > options.max_three_prime_complementarity
                || metrics.hairpin.max_three_prime_run > options.max_three_prime_complementarity
            {
                continue;
            }
            list.push(PrimerCandidate {
                primer,
                location,
                analysis: metrics,
            });
        }
    }
    let mut report=DesignReport {template_strand:molecule.template_strand,source_revision:source_revision.into(),options:options.clone(),pairs:vec![],evaluated_pairs:0,rejected_pairs:0,warnings:vec!["Anchored endpoint selection ranks the requested length interval under explicit heuristic filters. It does not establish binding uniqueness; inspect binding hits and independently simulate every PCR alternative before using a pair.".into(),"Ranking is a reproducible heuristic: summed absolute distance from target Tm plus pair Tm difference; no kinetics or empirical success probability is assigned.".into()]};
    let pair_work = forward.iter().try_fold(0usize, |sum, f| {
        reverse.iter().try_fold(sum, |sum, r| {
            sum.checked_add(
                f.primer
                    .ordered_sequence()
                    .len()
                    .saturating_mul(r.primer.ordered_sequence().len())
                    .saturating_mul(2),
            )
        })
    });
    if pair_work.is_none_or(|w| w > 20_000_000) {
        return Err(Error::Unsupported("primer-pair complementarity exceeds 20 million comparisons; narrow the length interval or shorten tails".into()));
    }
    for f in &forward {
        for r in &reverse {
            report.evaluated_pairs += 1;
            if f.primer.annealing.len() + r.primer.annealing.len() > target_len {
                report.rejected_pairs += 1;
                continue;
            }
            let tm_difference_c =
                (f.analysis.tm.as_ref().unwrap().tm_c - r.analysis.tm.as_ref().unwrap().tm_c).abs();
            let heterodimer = heterodimer(&f.primer, &r.primer)?;
            if tm_difference_c > options.max_tm_difference_c
                || heterodimer.max_three_prime_run > options.max_three_prime_complementarity
            {
                report.rejected_pairs += 1;
                continue;
            }
            let score = (f.analysis.tm.as_ref().unwrap().tm_c - options.target_tm_c).abs()
                + (r.analysis.tm.as_ref().unwrap().tm_c - options.target_tm_c).abs()
                + tm_difference_c;
            let pair = PrimerPair {
                forward: f.clone(),
                reverse: r.clone(),
                heterodimer,
                tm_difference_c,
                score,
                product_length: target_len + f.primer.tail.len() + r.primer.tail.len(),
            };
            let insert = report.pairs.partition_point(|p| p.score <= score);
            if insert < options.max_pairs {
                report.pairs.insert(insert, pair);
                report.pairs.truncate(options.max_pairs);
            }
        }
    }
    if report.pairs.is_empty() {
        report.warnings.push("No pair passed these endpoint, length, GC, Tm and contiguous-complementarity filters. Relax criteria deliberately or choose another target.".into());
    }
    Ok(report)
}
