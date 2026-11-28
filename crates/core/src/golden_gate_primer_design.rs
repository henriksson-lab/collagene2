//! Explicit source-interval Type IIS PCR preparation proposals. No target-directed DNA.
pub use crate::golden_gate::GoldenGateEnzyme;
use crate::{
    analysis::{self, Primer},
    assembly_recipe::{self, FragmentSpec, Preparation},
    catalog_digest::{self, CatalogSelection},
    chemistry::{self, Duplex, End, Polarity},
    verification::PlanInput,
};
use collagene_model::{Alphabet, Error, Location, Molecule, Sequence, Topology};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
pub const ENGINE_VERSION: &str = "golden-gate-primer-design-1";
const MAX_ROWS: usize = 32;
const MAX_BASES: usize = 8 << 20;
const MAX_WORK: usize = 100_000_000;
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Fragment{pub role:String,pub input_id:String,pub start:usize,pub end:usize,pub reverse:bool,pub forward_length:usize,pub reverse_length:usize,pub right_overhang:String}}
dto! {pub struct Request{pub name:String,pub fragments:Vec<Fragment>,pub enzyme:GoldenGateEnzyme,pub clamp:String,pub spacer:String}}
dto! {pub struct TailSegments{pub clamp:String,pub motif:String,pub spacer:String,pub fusion:String}}
impl TailSegments {
    pub fn sequence(&self) -> String {
        format!("{}{}{}{}", self.clamp, self.motif, self.spacer, self.fusion)
    }
}
dto! {pub struct Row{pub role:String,pub input_id:String,pub input_revision:String,pub forward:Primer,pub reverse:Primer,pub forward_tail:TailSegments,pub reverse_tail:TailSegments,pub selection:CatalogSelection,pub fragment_index:usize,pub oriented_top:String,pub left:End,pub right:End}}
dto! {pub struct Report{pub context_revision:String,pub source_revisions:Vec<String>,pub request:Request,pub fragments:Vec<FragmentSpec>,pub rows:Vec<Row>,pub warnings:Vec<String>}}
fn bad(s: &str) -> Error {
    Error::Unsupported(format!("Golden Gate primer planning: {s}"))
}
fn text(s: &str, n: usize) -> Result<(), Error> {
    if s.len() > n || s.trim().is_empty() || s.chars().any(char::is_control) {
        Err(bad("invalid bounded identifier"))
    } else {
        Ok(())
    }
}
fn dna(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b"ACGT".contains(&b))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn rc(s: &str) -> Result<String, Error> {
    Ok(Sequence::dna(s)?.reverse_complement()?.as_str().into())
}
fn dimensions(r: &Request) -> Result<(usize, usize), Error> {
    text(&r.name, 1024)?;
    let p = r.enzyme.profile();
    let spacer = usize::try_from(p.top_cut)
        .map_err(|_| Error::TooLarge)?
        .checked_sub(p.recognition.len())
        .ok_or_else(|| bad("unsupported cut geometry"))?;
    let k = usize::try_from(p.bottom_cut - p.top_cut).map_err(|_| Error::TooLarge)?;
    if r.fragments.is_empty()
        || r.fragments.len() > MAX_ROWS
        || r.clamp.len() > 12
        || !dna(&r.clamp)
        || r.spacer.len() != spacer
        || !dna(&r.spacer)
    {
        return Err(bad(
            "requires1..32 rows, concrete uppercase clamp1..12 and exact profile-specific spacer",
        ));
    }
    let mut roles = BTreeSet::new();
    for f in &r.fragments {
        text(&f.role, 128)?;
        text(&f.input_id, 128)?;
        if !roles.insert(&f.role)
            || !(8..=80).contains(&f.forward_length)
            || !(8..=80).contains(&f.reverse_length)
            || f.right_overhang.len() != k
            || !dna(&f.right_overhang)
        {
            return Err(bad(
                "distinct roles, annealing lengths8..80 and concrete profile-length right fusion required",
            ));
        }
    }
    Ok((spacer, k))
}
pub(crate) fn physical_source(m: &Molecule) -> Result<(), Error> {
    if chemistry::has_unmapped_physical_evidence(m)
        || m.metadata.keys().any(|k| {
            k == "duplex_model"
                || k == "tracked_modifications"
                || k == "source:declared_termini"
                || k == "source:legacy:chemistry_unmapped"
                || k == "source:snapgene:chemistry_unmapped"
                || k.starts_with("assembly:intent:")
                || k.starts_with("topo:")
                || k.starts_with("variant:")
                || k.contains("recombination:")
                || k.contains("restriction_inspection:")
                || k.contains("restriction_catalog:")
        })
    {
        return Err(bad(
            "use original physical DNA sources, not detached reference/physical-intermediate claims",
        ));
    }
    if m.original_genbank.as_ref().is_some_and(|raw| {
        raw.lines()
            .filter_map(|l| {
                l.strip_prefix("COMMENT     ")
                    .or_else(|| l.strip_prefix("CC   "))
            })
            .any(|l| {
                matches!(
                    l.trim_end(),
                    "Collagene physical product v1: unverified sequence-only evidence"
                        | "Collagene reference-only v1: unverified sequence-only evidence"
                )
            })
    }) {
        return Err(bad(
            "imported physical/reference-only precursor is not an original PCR substrate",
        ));
    }
    Ok(())
}
/// Borrowed shape/aggregate-work guard before cloning, hashing or scientific scans.
pub fn source_preflight(
    inputs: &[(&str, &Molecule)],
    target: &Molecule,
    r: &Request,
) -> Result<(), Error> {
    let (_, k) = dimensions(r)?;
    if inputs.is_empty() || inputs.len() > MAX_ROWS {
        return Err(Error::TooLarge);
    }
    let mut models = inputs.iter().map(|(_, m)| *m).collect::<Vec<_>>();
    models.push(target);
    collagene_model::bounds::model_preflight(&models)?;
    let mut lookup = BTreeMap::new();
    let mut bases = target.sequence.len();
    if target.sequence.alphabet() != Alphabet::Dna {
        return Err(bad("independent target must be DNA"));
    }
    for &(id, m) in inputs {
        text(id, 128)?;
        if lookup.insert(id, m).is_some() {
            return Err(bad("duplicate input ID"));
        }
        bases = bases
            .checked_add(m.sequence.len())
            .filter(|n| *n <= MAX_BASES)
            .ok_or(Error::TooLarge)?;
    }
    let mut work = 0usize;
    let mut output = 0usize;
    for f in &r.fragments {
        let m = *lookup
            .get(f.input_id.as_str())
            .ok_or_else(|| bad("row has no bound original input"))?;
        let n = m.sequence.len();
        let l = f
            .end
            .checked_sub(f.start)
            .filter(|l| *l <= n && *l > k)
            .ok_or_else(|| {
                bad("interval must contain at most one source turn and exceed fusion length")
            })?;
        let annealed = l - if f.reverse { k } else { 0 };
        if f.start >= n
            || f.end > f.start.checked_add(n).ok_or(Error::TooLarge)?
            || m.topology == Topology::Linear && f.end > n
            || f.forward_length + f.reverse_length > annealed
        {
            return Err(bad(
                "invalid source interval or overlapping annealing domains",
            ));
        }
        if m.sequence.alphabet() != Alphabet::Dna {
            return Err(bad("PCR requires original DNA"));
        }
        physical_source(m)?;
        work = work
            .checked_add(
                n.checked_mul((f.forward_length + f.reverse_length + 16) * 8)
                    .ok_or(Error::TooLarge)?,
            )
            .filter(|n| *n <= MAX_WORK)
            .ok_or(Error::TooLarge)?;
        output = output
            .checked_add(l + 128)
            .filter(|n| *n <= 2 << 20)
            .ok_or(Error::TooLarge)?;
    }
    // Validation/linear scans occur only after aggregate shape/work reservation.
    for m in models {
        m.validate()?;
    }
    for f in &r.fragments {
        let m = lookup[f.input_id.as_str()];
        if !dna(m.sequence.as_str()) {
            return Err(bad(
                "ambiguous source DNA prevents complete unique-PCR/cleavage evidence",
            ));
        }
    }
    Ok(())
}
pub fn preview_revision(
    inputs: &[PlanInput],
    target: &Molecule,
    r: &Request,
) -> Result<String, Error> {
    if inputs.is_empty() || inputs.len() > MAX_ROWS {
        return Err(Error::TooLarge);
    }
    source_preflight(
        &inputs
            .iter()
            .map(|i| (i.id.as_str(), &i.molecule))
            .collect::<Vec<_>>(),
        target,
        r,
    )?;
    let t = assembly_recipe::bound_revisions(inputs, target)?;
    assembly_recipe::hash(&(
        ENGINE_VERSION,
        inputs
            .iter()
            .map(|i| (&i.id, &i.revision))
            .collect::<Vec<_>>(),
        t,
        r,
    ))
}
fn extract(m: &Molecule, start: usize, end: usize) -> Result<String, Error> {
    let n = m.sequence.len();
    let len = end - start;
    let first = len.min(n - start % n);
    let loc = if first == len {
        Location::range(start % n, start % n + len)
    } else {
        Location::Join(vec![
            Location::range(start % n, n),
            Location::range(0, len - first),
        ])
    };
    loc.extract(&m.sequence)
}
fn tail(r: &Request, fusion: String) -> TailSegments {
    TailSegments {
        clamp: r.clamp.clone(),
        motif: r.enzyme.profile().recognition,
        spacer: r.spacer.clone(),
        fusion,
    }
}
pub fn design(inputs: &[PlanInput], target: &Molecule, r: &Request) -> Result<Report, Error> {
    let context_revision = preview_revision(inputs, target, r)?;
    let catalog = catalog_digest::builtin_catalog()?;
    let profile = r.enzyme.profile();
    let selection = catalog_digest::select(&catalog, &[profile.name.to_ascii_lowercase()])?;
    let enzymes = selection.resolved_enzymes(&catalog)?;
    if enzymes.len() != 1 || assembly_recipe::hash(&enzymes[0])? != assembly_recipe::hash(&profile)?
    {
        return Err(bad(
            "embedded catalog geometry differs from declared enzyme profile",
        ));
    }
    let mut rows = Vec::new();
    let mut fragments = Vec::new();
    for f in &r.fragments {
        let input = inputs
            .iter()
            .find(|i| i.id == f.input_id)
            .ok_or_else(|| bad("missing original"))?;
        let (fragment, row) = design_row(input, r, f, "", &selection, &enzymes)?;
        fragments.push(fragment);
        rows.push(row);
    }
    let report=Report{context_revision,source_revisions:inputs.iter().map(|i|i.revision.clone()).collect(),request:r.clone(),fragments,rows,warnings:vec!["Source-bound exact-PCR/cut geometry only, not empirical amplification, cleavage efficiency, methylation clearance or ligation-mixture uniqueness. Clamp/spacer are operator choices, not kit recommendations.".into(),"Right fusion is specified in the desired assembly-top reference direction; the actual protruding right-bottom strand is its reverse complement. Reverse rows use the original bottom-strand geometry, not reverse-complementing a full selected top fragment after digestion.".into(),"The independent target is bound but never determines oligo sequence. Freshly verify the returned preparations and full ordered ligation against that target; existing-site rows may remain separate.".into()]};
    report.validate()?;
    Ok(report)
}
/// Shared actual original-source PCR and physical digest reconstruction. The
/// old public design path passes an empty insert and keeps its serialized DTOs.
pub(crate) fn design_row(
    input: &PlanInput,
    r: &Request,
    f: &Fragment,
    inserted: &str,
    selection: &CatalogSelection,
    enzymes: &[analysis::Enzyme],
) -> Result<(FragmentSpec, Row), Error> {
    design_row_with_flanks(input, r, f, "", inserted, selection, enzymes)
}
/// Same actual PCR/digest materializer with explicitly introduced assembly-top
/// prefix/suffix. Old public paths pass an empty prefix and retain their DTOs.
pub(crate) fn design_row_with_flanks(
    input: &PlanInput,
    r: &Request,
    f: &Fragment,
    prefix: &str,
    inserted: &str,
    selection: &CatalogSelection,
    enzymes: &[analysis::Enzyme],
) -> Result<(FragmentSpec, Row), Error> {
    let (_, k) = dimensions(r)?;
    if prefix.len() > 200
        || !prefix.bytes().all(|b| b"ACGT".contains(&b))
        || inserted.len() > 200
        || !inserted.bytes().all(|b| b"ACGT".contains(&b))
    {
        return Err(bad("literal outgoing insert must be bounded concrete DNA"));
    }
    let profile = r.enzyme.profile();
    let m = &input.molecule;
    let core = extract(m, f.start, f.end)?;
    let annealed_end = f
        .end
        .checked_sub(if f.reverse { k } else { 0 })
        .ok_or_else(|| bad("source interval is shorter than its oriented fusion"))?;
    let annealed = extract(m, f.start, annealed_end)?;
    if f.forward_length
        .checked_add(f.reverse_length)
        .is_none_or(|n| n > annealed.len())
    {
        return Err(bad(
            "introduced flanks cannot replace overlapping or missing original annealing domains",
        ));
    }
    let mut ft = tail(
        r,
        if f.reverse {
            rc(&f.right_overhang)?
        } else {
            String::new()
        },
    );
    let mut rt = tail(
        r,
        if f.reverse {
            rc(&core[core.len() - k..])?
        } else {
            rc(&f.right_overhang)?
        },
    );
    if f.reverse {
        rt.fusion = format!("{prefix}{}", rt.fusion);
    } else {
        ft.fusion.push_str(prefix);
    }
    if !inserted.is_empty() {
        if f.reverse {
            ft.fusion.push_str(&rc(inserted)?);
        } else {
            rt.fusion.push_str(&rc(inserted)?);
        }
    }
    let forward = Primer::new(
        format!("{}-forward", f.role),
        &annealed[..f.forward_length],
        &ft.sequence(),
    )?;
    let reverse = Primer::new(
        format!("{}-reverse", f.role),
        &rc(&annealed[annealed.len() - f.reverse_length..])?,
        &rt.sequence(),
    )?;
    if forward.annealing.len() + forward.tail.len() > 200
        || reverse.annealing.len() + reverse.tail.len() > 200
    {
        return Err(bad(
            "complete ordered oligo including literal insert exceeds200 bases",
        ));
    }
    // Reject alternate binding sites before the PCR materializer can enumerate/copy them.
    let fwd = analysis::find_matches(m, forward.annealing.as_str(), false)?;
    let frev = analysis::find_matches(m, &rc(forward.annealing.as_str())?, false)?;
    let rplus = analysis::find_matches(m, reverse.annealing.as_str(), false)?;
    let rminus = analysis::find_matches(m, &rc(reverse.annealing.as_str())?, false)?;
    if fwd.len() != 1
        || !frev.is_empty()
        || !rplus.is_empty()
        || rminus.len() != 1
        || fwd[0].start != f.start % m.sequence.len()
        || rminus[0].start != (annealed_end - f.reverse_length) % m.sequence.len()
    {
        return Err(bad(
            "anchored primers do not have a unique exact original-source geometry",
        ));
    }
    let products = analysis::simulate_pcr_geometry(m, &forward, &reverse)?;
    if products.len() != 1 {
        return Err(bad(
            "exact PCR does not uniquely reconstruct the declared interval",
        ));
    }
    let p = &products[0].product;
    let expected = format!("{}{}{}", forward.tail, annealed, rc(&reverse.tail)?);
    if p.sequence.as_str() != expected {
        return Err(bad(
            "fresh PCR differs from independently declared primer/source geometry",
        ));
    }
    let amplicon = Molecule::dna("proposed PCR", &expected, Topology::Linear)?;
    let sites = analysis::restriction_sites(&amplicon, &profile)?;
    let left = r.clamp.len() + profile.recognition.len() + r.spacer.len();
    let right = left + prefix.len() + core.len() + inserted.len();
    if sites.len() != 2
        || sites.iter().any(|s| !s.certain)
        || !sites
            .iter()
            .any(|s| !s.reverse && s.top_cut == left && s.bottom_cut == left + k)
        || !sites
            .iter()
            .any(|s| s.reverse && s.top_cut == right && s.bottom_cut == right + k)
    {
        return Err(bad(
            "requires exactly the two eligible added Type IIS cuts; internal/extra sites are unsupported",
        ));
    }
    let cut = chemistry::digest(&Duplex::exact_pcr(&p.sequence)?, enzymes)?;
    let selected = cut
        .iter()
        .enumerate()
        .filter(|(_, d)| {
            d.top_start == left && d.duplex.top.len() == prefix.len() + core.len() + inserted.len()
        })
        .collect::<Vec<_>>();
    if selected.len() != 1 {
        return Err(bad("digest does not isolate exactly one intended core"));
    }
    let (index, fragment) = selected[0];
    let oriented = if f.reverse {
        fragment.duplex.reverse()?
    } else {
        fragment.duplex.clone()
    };
    let wanted = format!(
        "{}{}{}",
        prefix,
        if f.reverse { rc(&core)? } else { core.clone() },
        inserted
    );
    let (le, re) = oriented.ends()?;
    if oriented.top.as_str() != wanted
        || le.polarity != Polarity::FivePrime
        || re.polarity != Polarity::FivePrime
        || le.overhang != wanted[..k]
        || re.overhang != rc(&f.right_overhang)?
    {
        return Err(bad(
            "actual staggered duplex/orientation/fusion differs from declared source core",
        ));
    }
    let fragment = FragmentSpec {
        generation: None,
        role: f.role.clone(),
        input_id: f.input_id.clone(),
        reverse: f.reverse,
        preparation: Preparation::PcrCatalogDigest {
            forward: forward.clone(),
            reverse: reverse.clone(),
            selection: selection.clone(),
            fragment_index: index,
            repair: None,
            template_strand: crate::chemistry::Strand::Top,
        },
    };
    let row = Row {
        role: f.role.clone(),
        input_id: f.input_id.clone(),
        input_revision: input.revision.clone(),
        forward,
        reverse,
        forward_tail: ft,
        reverse_tail: rt,
        selection: selection.clone(),
        fragment_index: index,
        oriented_top: wanted,
        left: le,
        right: re,
    };
    Ok((fragment, row))
}
impl Report {
    /// Bounded readonly consistency. Imported reports never authorize physical preparation.
    pub fn validate(&self) -> Result<(), Error> {
        let (_, k) = dimensions(&self.request)?;
        let n = self.request.fragments.len();
        if !hex(&self.context_revision)
            || self.rows.len() != n
            || self.fragments.len() != n
            || self.source_revisions.is_empty()
            || self.source_revisions.len() > MAX_ROWS
            || self.source_revisions.iter().any(|s| !hex(s))
            || self.warnings.len() > 32
            || self.warnings.iter().any(|s| s.len() > 8192)
        {
            return Err(bad("invalid report shape"));
        }
        let mut bases = 0usize;
        for ((f, row), spec) in self
            .request
            .fragments
            .iter()
            .zip(&self.rows)
            .zip(&self.fragments)
        {
            let len = f
                .end
                .checked_sub(f.start)
                .filter(|n| *n > k && *n <= 2 << 20)
                .ok_or(Error::TooLarge)?;
            bases = bases
                .checked_add(row.oriented_top.len())
                .filter(|n| *n <= 2 << 20)
                .ok_or(Error::TooLarge)?;
            if row.oriented_top.len() != len
                || !dna(&row.oriented_top)
                || row.role != f.role
                || row.input_id != f.input_id
                || !hex(&row.input_revision)
                || !self.source_revisions.contains(&row.input_revision)
                || spec.role != f.role
                || spec.input_id != f.input_id
                || spec.reverse != f.reverse
                || row.fragment_index != 1
            {
                return Err(bad("row identity/geometry changed"));
            }
            row.selection.validate()?;
            let profile = self.request.enzyme.profile();
            if row.selection.enzymes.len() != 1 {
                return Err(bad("exact single-enzyme catalog pin required"));
            }
            let selected = &row.selection.enzymes[0];
            if selected.id != profile.name.to_ascii_lowercase()
                || selected.enzyme.name != profile.name
                || selected.enzyme.recognition != profile.recognition
                || selected.enzyme.top_cut != profile.top_cut
                || selected.enzyme.bottom_cut != profile.bottom_cut
            {
                return Err(bad("selected catalog enzyme geometry changed"));
            }
            for (t, p, length, suffix) in [
                (&row.forward_tail, &row.forward, f.forward_length, "forward"),
                (&row.reverse_tail, &row.reverse, f.reverse_length, "reverse"),
            ] {
                if p.name != format!("{}-{suffix}", f.role)
                    || p.annealing.len() != length
                    || p.annealing.alphabet() != Alphabet::Dna
                    || !dna(p.annealing.as_str())
                    || p.tail.len() > 64
                    || t.clamp != self.request.clamp
                    || t.motif != self.request.enzyme.profile().recognition
                    || t.spacer != self.request.spacer
                    || t.fusion.len() > k
                    || p.tail != t.sequence()
                {
                    return Err(bad("oligo/tail segmentation changed"));
                }
            }
            let forward_fusion = if f.reverse {
                rc(&f.right_overhang)?
            } else {
                String::new()
            };
            let reverse_fusion = if f.reverse {
                row.oriented_top[..k].to_string()
            } else {
                rc(&f.right_overhang)?
            };
            if row.forward_tail.fusion != forward_fusion
                || row.reverse_tail.fusion != reverse_fusion
                || row.left.chemistry != chemistry::Terminus::CLEAVED
                || row.right.chemistry != chemistry::Terminus::CLEAVED
                || row.left.polarity != Polarity::FivePrime
                || row.right.polarity != Polarity::FivePrime
                || row.left.overhang != row.oriented_top[..k]
                || row.right.overhang != rc(&f.right_overhang)?
            {
                return Err(bad("oriented fusion evidence changed"));
            }
            let Preparation::PcrCatalogDigest {
                forward,
                reverse,
                selection,
                fragment_index,
                repair: None,
                template_strand: crate::chemistry::Strand::Top,
                ..
            } = &spec.preparation
            else {
                return Err(bad(
                    "requires exact PCR then pinned catalog digest, no repair",
                ));
            };
            if *fragment_index != row.fragment_index
                || assembly_recipe::hash(forward)? != assembly_recipe::hash(&row.forward)?
                || assembly_recipe::hash(reverse)? != assembly_recipe::hash(&row.reverse)?
                || assembly_recipe::hash(selection)? != assembly_recipe::hash(&row.selection)?
            {
                return Err(bad("preparation disagrees with row"));
            }
        }
        assembly_recipe::hash(self)?;
        Ok(())
    }
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        self.validate()?;
        if self.context_revision != expected {
            return Err(bad("source/target/request context changed"));
        }
        Ok(())
    }
}
