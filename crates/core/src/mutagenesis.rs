//! Independently implemented, explicit primer-encoded substitutions and inverse-PCR plans.
//! Sequence compatibility is conditional evidence; no mismatch kinetics or kit outcome is inferred.
use crate::{
    Alphabet, Error, Location, Molecule, Sequence, Topology,
    analysis::{self, Primer},
    chemistry,
    model::MAX_SEQUENCE_LEN,
    primers::{self, BindingHit, BindingOptions, BindingReport},
    run_material::MaterialOrigin,
    verification::{
        self, DesignPlan, Operation, OrientedSource, PlanInput, Source, VerificationReport,
    },
};
use serde::{Deserialize, Serialize};
use std::io;

const MAX_BASES: usize = 32_000_000;
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct PcrOptions {
    pub binding: BindingOptions,
    pub max_pair_checks: usize,
    pub max_products: usize,
    pub max_product_bases: usize,
}
impl Default for PcrOptions {
    fn default() -> Self {
        Self {
            binding: BindingOptions {
                max_mismatches: 1,
                min_three_prime_exact: 1,
                ..Default::default()
            },
            max_pair_checks: 10_000,
            max_products: 1000,
            max_product_bases: MAX_BASES,
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IntroducedChange {
    pub source_position: usize,
    pub product_position: usize,
    pub template_base: char,
    pub product_base: char,
    pub template_ambiguous: bool,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MutantPcrProduct {
    pub sequence: Sequence,
    pub forward_start: usize,
    pub reverse_start: usize,
    pub template_length: usize,
    pub crosses_origin: bool,
    pub forward_primer_index: usize,
    pub reverse_primer_index: usize,
    pub forward_tail_length: usize,
    pub reverse_tail_length: usize,
    pub forward_binding: BindingHit,
    pub reverse_binding: BindingHit,
    pub substitutions: Vec<IntroducedChange>,
    pub overlapping_domains: bool,
    /// Concrete product bases and concrete binding evidence, not a claim of amplification.
    pub certain: bool,
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PcrReport {
    pub source_revision: String,
    pub options: PcrOptions,
    pub bindings: [BindingReport; 2],
    pub products: Vec<MutantPcrProduct>,
    pub pair_checks: usize,
    /// All placements and candidate pairs were visited within the configured budgets.
    pub search_complete: bool,
    /// Also requires no unsupported overlap/heteroduplex geometry or unresolved base evidence.
    pub enumeration_complete: bool,
    pub unsupported_geometries: usize,
    pub conflicting_overlaps: usize,
    pub ambiguous_bindings: usize,
    pub ambiguous_products: usize,
    pub rejected_geometries: usize,
    pub warnings: Vec<String>,
}
fn unsupported(message: &str) -> Error {
    Error::Unsupported(message.into())
}
fn checked_source(source: &Molecule) -> Result<(), Error> {
    source.validate()?;
    if source.sequence.alphabet() != Alphabet::Dna {
        return Err(unsupported("mutagenesis requires a DNA reference"));
    }
    if chemistry::has_unmapped_physical_evidence(source) {
        return Err(unsupported(
            "PCR on imported physical or declared non-duplex evidence requires its originating chemistry graph",
        ));
    }
    Ok(())
}
fn checked_primer(primer: &Primer) -> Result<Primer, Error> {
    if primer.annealing.alphabet() != Alphabet::Dna {
        return Err(unsupported("mutagenic oligos require DNA"));
    }
    let normalized = Primer::new(&primer.name, primer.annealing.as_str(), &primer.tail)?;
    if normalized.name.len() > 1024
        || normalized.name.chars().any(char::is_control)
        || normalized.ordered_sequence().len() > 200
    {
        return Err(unsupported(
            "oligos require a bounded name and at most 200 ordered DNA bases",
        ));
    }
    if normalized
        .ordered_sequence()
        .bytes()
        .any(|b| analysis::base_mask(b).count_ones() != 1)
    {
        return Err(unsupported(
            "degenerate oligo mixtures cannot establish a concrete mutagenic PCR product",
        ));
    }
    Ok(normalized)
}
fn rc(text: &str) -> Result<String, Error> {
    if text.is_empty() {
        Ok(String::new())
    } else {
        Ok(Sequence::dna(text)?.reverse_complement()?.as_str().into())
    }
}

/// Enumerate both supplied oligos in both physical roles, including same-oligo pairs.
/// Mismatch positions are recorded in each oligo's own 5′→3′ coordinates by the binding module.
/// For non-overlapping domains, final primer-derived bases replace their template counterparts.
/// Consistent overlaps are previewed but remain unsupported for complete PCR verification;
/// incompatible overlaps do not acquire an invented repaired duplex sequence.
pub fn simulate_pcr(
    source: &Molecule,
    forward: &Primer,
    reverse: &Primer,
    options: &PcrOptions,
) -> Result<PcrReport, Error> {
    checked_source(source)?;
    simulate_pcr_checked(source, forward, reverse, options)
}
/// Explicit graph-state entry point. Molecule metadata cannot establish this state;
/// the caller must resolve a previously reconstructed duplex from the operation graph.
pub fn simulate_pcr_reconstructed(
    source: &Molecule,
    substrate: &chemistry::Duplex,
    forward: &Primer,
    reverse: &Primer,
    options: &PcrOptions,
) -> Result<PcrReport, Error> {
    source.validate()?;
    substrate.validate()?;
    if source.sequence.alphabet() != Alphabet::Dna
        || substrate.top != source.sequence
        || substrate.topology != source.topology
        || substrate.bottom_start != 0
        || substrate.bottom.len() != substrate.top.len()
        || !substrate.nicks.is_empty()
    {
        return Err(unsupported(
            "mismatched PCR needs matching, intact, fully paired reconstructed template strands",
        ));
    }
    simulate_pcr_checked(source, forward, reverse, options)
}
fn simulate_pcr_checked(
    source: &Molecule,
    forward: &Primer,
    reverse: &Primer,
    options: &PcrOptions,
) -> Result<PcrReport, Error> {
    if options.max_pair_checks == 0
        || options.max_pair_checks > 1_000_000
        || options.max_products == 0
        || options.max_products > 10_000
        || options.max_product_bases == 0
        || options.max_product_bases > MAX_BASES
    {
        return Err(unsupported(
            "invalid mutagenic PCR pair, product or materialization budget",
        ));
    }
    let oligos = [checked_primer(forward)?, checked_primer(reverse)?];
    let bindings = [
        primers::find_binding_sites(source, &oligos[0], &options.binding)?,
        primers::find_binding_sites(source, &oligos[1], &options.binding)?,
    ];
    let search_complete = !bindings.iter().any(|b| b.truncated);
    let ambiguous_bindings = bindings
        .iter()
        .flat_map(|b| &b.hits)
        .filter(|h| !h.ambiguous_positions.is_empty())
        .count();
    let mut report = PcrReport { source_revision: source.revision(), options: options.clone(), bindings,
        products: vec![], pair_checks: 0, search_complete, enumeration_complete: false,
        unsupported_geometries: 0, conflicting_overlaps: 0, ambiguous_bindings, ambiguous_products: 0,
        rejected_geometries: 0, warnings: vec![
            "Primer-domain mismatches are a declared compatibility rule, not a probability of extension. No mismatch thermodynamics, polymerase kinetics or experimental outcome is evaluated.".into(),
            "Both oligos are enumerated in both roles, including same-oligo products; names do not constrain orientation. Searches cover this reference only and at most one template turn.".into(),
            "Synthetic tails and actual annealing-domain oligo bases enter the final product. Consistent overlapping domains are sequence previews only; conflicting overlaps/heteroduplexes and overlapping whole-plasmid mutagenesis chemistry are not repaired or silently modeled.".into(),
        ] };
    let identical_oligos =
        oligos[0].annealing == oligos[1].annealing && oligos[0].tail == oligos[1].tail;
    let mut total_bases = 0usize;
    let n = source.sequence.len();
    let input = source.sequence.as_str().as_bytes();
    'pairs: for (fi, f) in oligos.iter().enumerate() {
        if identical_oligos && fi == 1 {
            continue;
        }
        for (ri, r) in oligos.iter().enumerate() {
            if identical_oligos && ri == 1 {
                continue;
            }
            let reverse_domain = rc(r.annealing.as_str())?;
            let reverse_tail = rc(&r.tail)?;
            for fh in report.bindings[fi].hits.iter().filter(|h| !h.reverse) {
                for rh in report.bindings[ri].hits.iter().filter(|h| h.reverse) {
                    if report.pair_checks == options.max_pair_checks {
                        report.search_complete = false;
                        report.warnings.push("Candidate-pair checking reached its budget; additional products have not been excluded.".into());
                        break 'pairs;
                    }
                    report.pair_checks += 1;
                    let end = rh.start + rh.length;
                    let length = if source.topology == Topology::Circular {
                        let d = (end + n - fh.start) % n;
                        if d == 0 { n } else { d }
                    } else if end > fh.start {
                        end - fh.start
                    } else {
                        report.rejected_geometries += 1;
                        continue;
                    };
                    if length < fh.length.max(rh.length) || length > n {
                        report.rejected_geometries += 1;
                        continue;
                    }
                    let overlapping = length < fh.length + rh.length;
                    if overlapping {
                        report.unsupported_geometries += 1;
                    }
                    let product_length = f.tail.len() + length + reverse_tail.len();
                    if product_length > MAX_SEQUENCE_LEN {
                        report.unsupported_geometries += 1;
                        continue;
                    }
                    let right = length - rh.length;
                    if reverse_domain.bytes().enumerate().any(|(j, b)| {
                        right + j < fh.length && f.annealing.as_str().as_bytes()[right + j] != b
                    }) {
                        report.conflicting_overlaps += 1;
                        continue;
                    }
                    if report.products.len() == options.max_products
                        || total_bases.saturating_add(product_length) > options.max_product_bases
                    {
                        report.search_complete = false;
                        report.warnings.push("Product count/base materialization reached its budget; additional products have not been excluded.".into());
                        break 'pairs;
                    }
                    let mut template: Vec<u8> =
                        (0..length).map(|j| input[(fh.start + j) % n]).collect();
                    template[..fh.length].copy_from_slice(f.annealing.as_str().as_bytes());
                    template[right..].copy_from_slice(reverse_domain.as_bytes());
                    let substitutions = template
                        .iter()
                        .enumerate()
                        .filter_map(|(j, &base)| {
                            let p = (fh.start + j) % n;
                            (base != input[p]).then_some(IntroducedChange {
                                source_position: p,
                                product_position: f.tail.len() + j,
                                template_base: input[p] as char,
                                product_base: base as char,
                                template_ambiguous: analysis::base_mask(input[p]).count_ones() != 1,
                            })
                        })
                        .collect();
                    let mut bases = String::with_capacity(product_length);
                    bases.push_str(&f.tail);
                    bases.push_str(std::str::from_utf8(&template).expect("validated ASCII DNA"));
                    bases.push_str(&reverse_tail);
                    let certain = bases
                        .bytes()
                        .all(|b| analysis::base_mask(b).count_ones() == 1)
                        && fh.ambiguous_positions.is_empty()
                        && rh.ambiguous_positions.is_empty();
                    if !certain {
                        report.ambiguous_products += 1;
                    }
                    total_bases += product_length;
                    report.products.push(MutantPcrProduct {
                        sequence: Sequence::dna(&bases)?,
                        forward_start: fh.start,
                        reverse_start: rh.start,
                        template_length: length,
                        crosses_origin: fh.start + length > n,
                        forward_primer_index: fi,
                        reverse_primer_index: ri,
                        forward_tail_length: f.tail.len(),
                        reverse_tail_length: reverse_tail.len(),
                        forward_binding: fh.clone(),
                        reverse_binding: rh.clone(),
                        substitutions,
                        overlapping_domains: overlapping,
                        certain,
                    });
                }
            }
        }
    }
    report.enumeration_complete = report.search_complete
        && report.unsupported_geometries == 0
        && report.conflicting_overlaps == 0
        && report.ambiguous_bindings == 0
        && report.ambiguous_products == 0;
    if report.unsupported_geometries > 0 {
        report.warnings.push(format!(
            "{} unsupported overlap or document-size geometries prevent a complete product claim.",
            report.unsupported_geometries
        ));
    }
    if report.conflicting_overlaps > 0 {
        report.warnings.push(format!("{} overlapping primer pairs encode conflicting bases; no repaired product is asserted.",report.conflicting_overlaps));
    }
    if report.ambiguous_bindings > 0 || report.ambiguous_products > 0 {
        report.warnings.push("Ambiguous template evidence prevents concrete complete enumeration, even when a primer supplies a concrete final base.".into());
    }
    if !report.search_complete {
        report.warnings.push("Search is truncated. The displayed product set cannot establish uniqueness or absence of another product.".into());
    }
    Ok(report)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct InversePcrOptions {
    /// Half-open source interval to replace; end may pass the circular origin by at most one turn.
    pub edit_start: usize,
    pub edit_end: usize,
    pub inserted: String,
    pub forward_length: usize,
    pub reverse_length: usize,
    /// Initial inserted bases assigned to the reverse oligo as a reverse-complemented 5′ tail.
    pub reverse_tail_bases: usize,
}
impl Default for InversePcrOptions {
    fn default() -> Self {
        Self {
            edit_start: 0,
            edit_end: 0,
            inserted: String::new(),
            forward_length: 20,
            reverse_length: 20,
            reverse_tail_bases: 0,
        }
    }
}
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InversePcrDesign {
    pub context_revision: String,
    pub source_revision: String,
    pub target_revision: String,
    pub options: InversePcrOptions,
    pub forward: Primer,
    pub reverse: Primer,
    pub forward_location: Location,
    pub reverse_location: Location,
    pub retained_template_length: usize,
    pub plan: DesignPlan,
    pub verification: VerificationReport,
    pub warnings: Vec<String>,
}
/// Pure inverse-PCR plan and actual oligos, before any molecular verification.
/// The target is supplied independently and never copied from a PCR product.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub struct InversePcrPreparation {
    pub context_revision: String,
    pub source_revision: String,
    pub target_revision: String,
    pub options: InversePcrOptions,
    pub forward: Primer,
    pub reverse: Primer,
    pub forward_location: Location,
    pub reverse_location: Location,
    pub retained_template_length: usize,
    pub plan: DesignPlan,
    pub warnings: Vec<String>,
}
fn region(start: usize, count: usize, n: usize, reverse: bool) -> Location {
    let first = count.min(n - start);
    let mut parts = vec![Location::range(start, start + first)];
    if first < count {
        parts.push(Location::range(0, count - first));
    }
    let location = if parts.len() == 1 {
        parts.remove(0)
    } else {
        Location::Join(parts)
    };
    if reverse {
        Location::Complement(Box::new(location))
    } else {
        location
    }
}
struct SizeLimit(usize);
impl io::Write for SizeLimit {
    fn write(&mut self, data: &[u8]) -> io::Result<usize> {
        self.0 = self
            .0
            .checked_add(data.len())
            .filter(|n| *n <= 32 * 1024 * 1024)
            .ok_or_else(|| io::Error::other("mutagenesis result exceeds wire limit"))?;
        Ok(data.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}
/// Bounded documentary request preflight. Saved-origin chemistry is reconstructed
/// only by the ordinary verifier, never recovered from source metadata.
fn inverse_source_preflight(
    source: &Molecule,
    origin: Option<&MaterialOrigin>,
) -> Result<(), Error> {
    crate::model::bounds::model_preflight(&[source])?;
    crate::run_material::material_origins_preflight(origin.into_iter())?;
    source.validate()?;
    if source.sequence.alphabet() != Alphabet::Dna || source.topology != Topology::Circular {
        return Err(unsupported("inverse PCR requires a circular DNA template"));
    }
    if let Some(origin) = origin {
        origin.validate_documentary()?;
        if origin.event().revision != source.revision() {
            return Err(unsupported("inverse-PCR saved template revision changed"));
        }
    } else {
        checked_source(source)?;
    }
    Ok(())
}
pub(crate) struct InverseParts {
    pub(crate) forward: Primer,
    pub(crate) reverse: Primer,
    pub(crate) forward_location: Location,
    pub(crate) reverse_location: Location,
    pub(crate) retained: usize,
}
pub(crate) fn inverse_parts(
    sequence: &Sequence,
    options: &InversePcrOptions,
) -> Result<InverseParts, Error> {
    // Bound the raw literal before normalization/allocation, not just the resulting oligos.
    if options.inserted.len() > 1024 {
        return Err(Error::TooLarge);
    }
    let n = sequence.len();
    let deleted = options
        .edit_end
        .checked_sub(options.edit_start)
        .ok_or_else(|| Error::Coordinates("edit end precedes start".into()))?;
    let terminal_caret = options.edit_start == n && options.edit_end == n;
    if (options.edit_start >= n && !terminal_caret)
        || deleted >= n
        || options.forward_length < 8
        || options.reverse_length < 8
        || options.forward_length > 80
        || options.reverse_length > 80
    {
        return Err(Error::Coordinates("replace less than one source turn, with bounded 8–80 base non-overlapping primer domains".into()));
    }
    let retained = n - deleted;
    if options.forward_length + options.reverse_length > retained {
        return Err(unsupported(
            "inverse-PCR annealing domains overlap within the retained source region",
        ));
    }
    let inserted = if options.inserted.trim().is_empty() {
        String::new()
    } else {
        Sequence::dna(&options.inserted)?.as_str().to_string()
    };
    if inserted
        .bytes()
        .any(|b| analysis::base_mask(b).count_ones() != 1)
        || options.reverse_tail_bases > inserted.len()
    {
        return Err(unsupported(
            "replacement tails require concrete DNA and a valid split position",
        ));
    }
    let forward_start = options.edit_end % n;
    let reverse_start = (options.edit_start + n - options.reverse_length) % n;
    let forward_location = region(forward_start, options.forward_length, n, false);
    let reverse_location = region(reverse_start, options.reverse_length, n, true);
    let forward = checked_primer(&Primer::new(
        "Mutagenesis forward",
        &forward_location.extract(sequence)?,
        &inserted[options.reverse_tail_bases..],
    )?)?;
    let reverse = checked_primer(&Primer::new(
        "Mutagenesis reverse",
        &reverse_location.extract(sequence)?,
        &rc(&inserted[..options.reverse_tail_bases])?,
    )?)?;
    if retained
        .checked_add(inserted.len())
        .is_none_or(|len| len > MAX_SEQUENCE_LEN)
    {
        return Err(Error::TooLarge);
    }
    Ok(InverseParts {
        forward,
        reverse,
        forward_location,
        reverse_location,
        retained,
    })
}
const INVERSE_ENGINE_VERSION: &str = "inverse-pcr-preparation-2";
fn inverse_context(
    source: &Molecule,
    target_revision: &str,
    origin: Option<&MaterialOrigin>,
    options: &InversePcrOptions,
) -> Result<String, Error> {
    crate::assembly_recipe::hash(&(
        INVERSE_ENGINE_VERSION,
        verification::ENGINE_VERSION,
        source.revision(),
        target_revision,
        origin,
        options,
    ))
}
/// Full source/origin/independent-target/options binding. No replay or PCR search.
pub fn inverse_pcr_revision(
    source: &Molecule,
    target: &Molecule,
    origin: Option<&MaterialOrigin>,
    options: &InversePcrOptions,
) -> Result<String, Error> {
    crate::model::bounds::model_preflight(&[source, target])?;
    // Serialize borrowed inputs under the aggregate envelope bound before
    // source/origin hashing, cloning, normalization or molecular work.
    serde_json::to_writer(&mut SizeLimit(65536), &(source, target, origin, options))
        .map_err(|_| Error::TooLarge)?;
    inverse_source_preflight(source, origin)?;
    target.validate()?;
    if target.topology != Topology::Circular || target.sequence.alphabet() != Alphabet::Dna {
        return Err(unsupported(
            "inverse PCR requires an independently supplied circular DNA target",
        ));
    }
    inverse_parts(&source.sequence, options)?;
    inverse_context(source, &target.revision(), origin, options)
}
/// Tail-encoded replacement, explicit phosphorylation and intended circular ligation.
/// The supplied target never chooses a PCR geometry. A selected saved origin is
/// reconstructed natively before PCR; no physical-source marker is removed.
pub fn compile_inverse_pcr(
    source: &Molecule,
    target: &Molecule,
    origin: Option<&MaterialOrigin>,
    options: &InversePcrOptions,
) -> Result<InversePcrPreparation, Error> {
    let context_revision = inverse_pcr_revision(source, target, origin, options)?;
    let InverseParts {
        forward,
        reverse,
        forward_location,
        reverse_location,
        retained,
    } = inverse_parts(&source.sequence, options)?;
    let mut steps = Vec::new();
    if let Some(origin) = origin {
        steps.push(match origin {
            MaterialOrigin::Run(origin) => Operation::ReconstructRunProduct {
                template: Source::Input("template".into()),
                origin: origin.clone(),
            },
            MaterialOrigin::Derived(origin) => Operation::ReconstructMaterial {
                template: Source::Input("template".into()),
                origin: origin.clone(),
            },
        });
    }
    let pcr_step = steps.len();
    steps.push(Operation::Pcr {
        template: if origin.is_some() {
            Source::Step(0)
        } else {
            Source::Input("template".into())
        },
        forward: forward.clone(),
        reverse: reverse.clone(),
        template_strand: crate::chemistry::Strand::Top,
    });
    steps.push(Operation::TreatEnds {
        template: Source::Step(pcr_step),
        treatment: chemistry::EndTreatment::Phosphorylate,
    });
    steps.push(Operation::Ligate {
        fragments: vec![OrientedSource {
            source: Source::Step(pcr_step + 1),
            reverse: false,
        }],
        circular: true,
        allow_nicked: false,
        options: Default::default(),
    });
    let plan = DesignPlan {
        name: format!("Inverse-PCR replacement in {}", source.name),
        target_revision: target.revision(),
        inputs: vec![PlanInput {
            id: "template".into(),
            revision: source.revision(),
            molecule: source.clone(),
        }],
        steps,
    };
    crate::run_material::plan_preflight(&plan)?;
    let mut warnings=vec![
        "Tail-encoded back-to-back inverse-PCR design: the replaced reference interval is omitted and the declared insertion is supplied by primer tails. No altered reference bases are smuggled in as a claimed PCR product.".into(),
        "The returned plan explicitly phosphorylates ideal newly synthesized unmodified PCR termini before intended sealed blunt circular ligation. This does not model a commercial reagent mixture, template digestion/removal, transformation, repair or empirical success.".into(),
        "Inspect fresh verification and all exact-PCR alternatives. A matching target is conditional sequence/chemistry consistency, not experimental confirmation or isolation.".into(),
    ];
    if origin.is_some() {
        warnings.push("The exact saved origin is replayed by fresh verification. Its inherited incompleteness and assumptions remain; nicked or unsupported circular strand geometry is not silently repaired.".into());
    }
    let prepared = InversePcrPreparation {
        context_revision,
        source_revision: source.revision(),
        target_revision: target.revision(),
        options: options.clone(),
        forward,
        reverse,
        forward_location,
        reverse_location,
        retained_template_length: retained,
        plan,
        warnings,
    };
    serde_json::to_writer(&mut SizeLimit(0), &prepared).map_err(|_| Error::TooLarge)?;
    Ok(prepared)
}
/// Compile the one current source-aware graph and perform fresh native verification.
pub fn design_inverse_pcr(
    source: &Molecule,
    target: &Molecule,
    origin: Option<&MaterialOrigin>,
    options: &InversePcrOptions,
) -> Result<InversePcrDesign, Error> {
    let prepared = compile_inverse_pcr(source, target, origin, options)?;
    let verification = verification::verify(&prepared.plan, target);
    let report = InversePcrDesign {
        context_revision: prepared.context_revision,
        source_revision: prepared.source_revision,
        target_revision: prepared.target_revision,
        options: prepared.options,
        forward: prepared.forward,
        reverse: prepared.reverse,
        forward_location: prepared.forward_location,
        reverse_location: prepared.reverse_location,
        retained_template_length: prepared.retained_template_length,
        plan: prepared.plan,
        verification,
        warnings: prepared.warnings,
    };
    serde_json::to_writer(&mut SizeLimit(0), &report).map_err(|_| Error::TooLarge)?;
    Ok(report)
}
// This borrowed validator checks documentary correspondence only. It does not
// replay material, enumerate PCR geometries, or adopt a historical assessment.
fn validate_inverse_fields(
    expected: &str,
    context: &str,
    source_revision: &str,
    target_revision: &str,
    options: &InversePcrOptions,
    forward: &Primer,
    reverse: &Primer,
    forward_location: &Location,
    reverse_location: &Location,
    retained: usize,
    plan: &DesignPlan,
) -> Result<(), Error> {
    let bad = || unsupported("inverse-PCR preparation binding or declared graph changed");
    if expected != context
        || plan.inputs.len() != 1
        || plan.inputs[0].id != "template"
        || plan.target_revision != target_revision
        || plan.inputs[0].revision != source_revision
    {
        return Err(bad());
    }
    crate::run_material::plan_preflight(plan)?;
    crate::model::bounds::model_preflight(&[&plan.inputs[0].molecule])?;
    let source = &plan.inputs[0].molecule;
    if source.revision() != source_revision
        || plan.name != format!("Inverse-PCR replacement in {}", source.name)
    {
        return Err(bad());
    }
    let origin = match plan.steps.first() {
        Some(Operation::ReconstructRunProduct {
            template: Source::Input(id),
            origin,
        }) if id == "template" => Some(MaterialOrigin::Run(origin.clone())),
        Some(Operation::ReconstructMaterial {
            template: Source::Input(id),
            origin,
        }) if id == "template" => Some(MaterialOrigin::Derived(origin.clone())),
        Some(Operation::Pcr {
            template: Source::Input(id),
            template_strand: crate::chemistry::Strand::Top,
            ..
        }) if id == "template" => None,
        _ => return Err(bad()),
    };
    inverse_source_preflight(source, origin.as_ref())?;
    if inverse_context(source, target_revision, origin.as_ref(), options)? != context {
        return Err(bad());
    }
    let parts = inverse_parts(&source.sequence, options)?;
    if crate::assembly_recipe::hash(&(
        forward,
        reverse,
        forward_location,
        reverse_location,
        retained,
    ))? != crate::assembly_recipe::hash(&(
        &parts.forward,
        &parts.reverse,
        &parts.forward_location,
        &parts.reverse_location,
        parts.retained,
    ))? {
        return Err(bad());
    }
    let i = usize::from(origin.is_some());
    if plan.steps.len() != i + 3 {
        return Err(bad());
    }
    let expected_template = if i == 0 {
        Source::Input("template".into())
    } else {
        Source::Step(0)
    };
    match &plan.steps[i] {
        Operation::Pcr {
            template,
            forward: f,
            reverse: r,
            template_strand: crate::chemistry::Strand::Top,
            ..
        } if crate::assembly_recipe::hash(&(template, f, r))?
            == crate::assembly_recipe::hash(&(expected_template, forward, reverse))? =>
        {
            ()
        }
        _ => return Err(bad()),
    }
    if !matches!(&plan.steps[i+1],Operation::TreatEnds{template:Source::Step(j),treatment:chemistry::EndTreatment::Phosphorylate} if *j==i)
    {
        return Err(bad());
    }
    match &plan.steps[i + 2] {
        Operation::Ligate {
            fragments,
            circular: true,
            allow_nicked: false,
            options,
        } if fragments.len() == 1
            && !fragments[0].reverse
            && matches!(&fragments[0].source,Source::Step(j) if *j==i+1)
            && crate::assembly_recipe::hash(options)?
                == crate::assembly_recipe::hash(&crate::ligation_outcomes::Options::default())? =>
        {
            ()
        }
        _ => return Err(bad()),
    }
    crate::run_material::plan_preflight(plan)
}
impl InversePcrPreparation {
    /// Readonly receipt/shape check, not molecular verification or source attestation.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        serde_json::to_writer(&mut SizeLimit(0), self).map_err(|_| Error::TooLarge)?;
        validate_inverse_fields(
            expected,
            &self.context_revision,
            &self.source_revision,
            &self.target_revision,
            &self.options,
            &self.forward,
            &self.reverse,
            &self.forward_location,
            &self.reverse_location,
            self.retained_template_length,
            &self.plan,
        )
    }
}
impl InversePcrDesign {
    /// Readonly receipt/shape check. Fresh verification still determines material.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        serde_json::to_writer(&mut SizeLimit(0), self).map_err(|_| Error::TooLarge)?;
        validate_inverse_fields(
            expected,
            &self.context_revision,
            &self.source_revision,
            &self.target_revision,
            &self.options,
            &self.forward,
            &self.reverse,
            &self.forward_location,
            &self.reverse_location,
            self.retained_template_length,
            &self.plan,
        )?;
        if self.verification.plan_revision != verification::plan_revision(&self.plan)
            || self.verification.target_revision != self.target_revision
            || self.verification.engine_version != verification::ENGINE_VERSION
            || self.verification.products.len() > self.plan.steps.len()
        {
            return Err(unsupported(
                "inverse-PCR verification receipt differs from the declared plan",
            ));
        }
        Ok(())
    }
}
