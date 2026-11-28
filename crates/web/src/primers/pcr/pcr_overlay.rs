//! Transient projection of selected native PCR evidence or captured selections.
//! No template scanning, candidate search, or experimental success inference.
use crate::primers::primer_pair_model::PairDraft;
use collagene_core::{
    cloning::pcr::{mutagenesis::PcrReport, primers::BindingHit},
    material::lineage::EventBinding,
    sequence::analysis::{Primer, base_mask},
    verification::pcr::PcrDocument,
};
use collagene_model::{Alphabet, Location, Molecule, Sequence, Topology};
use std::sync::Arc;

pub const MAX_SPANS: usize = 1024;
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Kind {
    ExactProduct,
    MismatchProduct,
    CapturedSelection,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Role {
    Forward,
    Reverse,
    Amplified,
    Mismatch,
    Ambiguous,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Span {
    pub start: usize,
    pub end: usize,
    pub reverse: bool,
    pub role: Role,
    pub label: String,
}
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Overlay {
    pub source: EventBinding,
    pub epoch: u64,
    pub kind: Kind,
    pub spans: Vec<Span>,
    /// Human-readable synthetic 5′ tails, always separate from reference spans.
    pub tails: Vec<String>,
    pub warnings: Vec<String>,
    pub summary: String,
}
#[derive(Clone)]
pub struct Publication {
    pub owner: String,
    pub overlay: Option<Arc<Overlay>>,
}
impl Overlay {
    pub fn is_current(&self, source: &EventBinding, epoch: u64) -> bool {
        &self.source == source && self.epoch == epoch
    }
    pub fn validate_for_length(&self, length: usize) -> bool {
        length > 0
            && self.spans.len() <= MAX_SPANS
            && self.spans.iter().all(|s| {
                s.start < s.end
                    && s.end <= length
                    && s.label.len() <= 1024
                    && !(self.kind == Kind::CapturedSelection && s.role == Role::Amplified)
            })
            && self.tails.len() <= 8
            && self.tails.iter().all(|s| s.len() <= 1024)
            && self.warnings.len() <= 64
            && self.warnings.iter().all(|s| s.len() <= 4096)
            && self.summary.len() <= 1024
    }
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn source_guard(m: &Molecule, source: &EventBinding) -> Result<(), String> {
    if source.document_id.is_empty()
        || source.document_id.len() > 128
        || !source
            .document_id
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b"_-.:".contains(&b))
        || !source.event_id.strip_prefix("event-").is_some_and(hex)
        || !hex(&source.revision)
    {
        return Err("Overlay requires the exact source document/event/revision.".into());
    }
    collagene_model::bounds::model_preflight(&[m]).map_err(|e| e.to_string())?;
    m.validate().map_err(|e| e.to_string())?;
    if m.sequence.alphabet() != Alphabet::Dna
        || m.sequence.is_empty()
        || m.revision() != source.revision
    {
        return Err("PCR overlay source changed or is not DNA.".into());
    }
    Ok(())
}
fn primer_guard(p: &Primer) -> Result<(), String> {
    if p.name.len() > 1024
        || p.name.chars().any(char::is_control)
        || p.annealing.alphabet() != Alphabet::Dna
        || p.annealing.is_empty()
        || p.annealing
            .len()
            .checked_add(p.tail.len())
            .is_none_or(|n| n > 200)
    {
        return Err("Overlay primer exceeds the bounded 200-base oligo profile.".into());
    }
    if !p.tail.is_empty() && Sequence::dna(&p.tail).map_err(|e| e.to_string())?.as_str() != p.tail {
        return Err("Overlay tail must be canonical DNA.".into());
    }
    Ok(())
}
fn new(source: &EventBinding, epoch: u64, kind: Kind) -> Overlay {
    Overlay {
        source: source.clone(),
        epoch,
        kind,
        spans: vec![],
        tails: vec![],
        warnings: vec![],
        summary: String::new(),
    }
}
fn finish(overlay: Overlay, length: usize) -> Result<Overlay, String> {
    if !overlay.validate_for_length(length) {
        return Err("Overlay exceeds its coordinate/output bounds.".into());
    }
    Ok(overlay)
}
fn push(
    overlay: &mut Overlay,
    n: usize,
    start: usize,
    length: usize,
    reverse: bool,
    role: Role,
    label: &str,
) {
    let first = length.min(n - start);
    overlay.spans.push(Span {
        start,
        end: start + first,
        reverse,
        role,
        label: label.into(),
    });
    if first < length {
        overlay.spans.push(Span {
            start: 0,
            end: length - first,
            reverse,
            role,
            label: label.into(),
        });
    }
}
fn domain_guard(m: &Molecule, start: usize, length: usize) -> Result<(), String> {
    let n = m.sequence.len();
    if start >= n
        || length == 0
        || length > n
        || m.topology == Topology::Linear && start + length > n
    {
        return Err("PCR placement is outside the source reference.".into());
    }
    Ok(())
}
fn label(index: usize, p: &Primer) -> String {
    format!(
        "Oligo {} · {}",
        index + 1,
        p.name.chars().take(128).collect::<String>()
    )
}
fn tail(overlay: &mut Overlay, p: &Primer, index: usize, direction: &str) {
    if !p.tail.is_empty() {
        let text = format!(
            "Oligo {} · {direction} 5′ synthetic tail (not template): {}",
            index + 1,
            p.tail
        );
        if !overlay.tails.contains(&text) {
            overlay.tails.push(text);
        }
    }
}
fn rc(s: &str) -> Result<String, String> {
    if s.is_empty() {
        return Ok(String::new());
    }
    Sequence::dna(s)
        .and_then(|v| v.reverse_complement())
        .map(|v| v.as_str().into())
        .map_err(|e| e.to_string())
}
fn template_length(m: &Molecule, start: usize, end: usize) -> Option<usize> {
    let n = m.sequence.len();
    if m.topology == Topology::Circular {
        let distance = (end + n - start) % n;
        Some(if distance == 0 { n } else { distance })
    } else {
        end.checked_sub(start).filter(|&length| length > 0)
    }
}

pub fn exact(
    m: &Molecule,
    source: &EventBinding,
    epoch: u64,
    forward: &Primer,
    reverse: &Primer,
    product: &PcrDocument,
) -> Result<Overlay, String> {
    primer_guard(forward)?;
    primer_guard(reverse)?;
    collagene_model::bounds::model_preflight(&[&product.molecule]).map_err(|e| e.to_string())?;
    source_guard(m, source)?;
    let g = &product.geometry;
    if product.molecule.sequence != g.sequence
        || g.sequence.alphabet() != Alphabet::Dna
        || product.molecule.topology != Topology::Linear
        || product.findings.len() > 64
        || product.findings.iter().any(|s| s.len() > 4096)
        || product
            .molecule
            .metadata
            .get("pcr:template_revision")
            .is_some_and(|revision| revision != &source.revision)
    {
        return Err("Selected exact-PCR preview differs from its source/sequence.".into());
    }
    let primers = [forward, reverse];
    let n = m.sequence.len();
    let source_bases = m.sequence.as_str().as_bytes();
    let mut overlay = new(source, epoch, Kind::ExactProduct);
    let mut alternatives = 0;
    for (fi, f) in primers.iter().enumerate() {
        for (ri, r) in primers.iter().enumerate() {
            let fl = f.annealing.len();
            let rl = r.annealing.len();
            if domain_guard(m, g.forward_start, fl).is_err()
                || domain_guard(m, g.reverse_start, rl).is_err()
            {
                continue;
            }
            let Some(length) = template_length(m, g.forward_start, g.reverse_start + rl) else {
                continue;
            };
            if length > n || length < fl + rl || g.crosses_origin != (g.forward_start + length > n)
            {
                continue;
            }
            let rd = rc(r.annealing.as_str())?;
            let rt = rc(&r.tail)?;
            if !f.annealing.as_str().bytes().enumerate().all(|(i, b)| {
                b == source_bases[(g.forward_start + i) % n] && base_mask(b).count_ones() == 1
            }) || !rd.bytes().enumerate().all(|(i, b)| {
                b == source_bases[(g.reverse_start + i) % n] && base_mask(b).count_ones() == 1
            }) {
                continue;
            }
            let sequence = g.sequence.as_str().as_bytes();
            if sequence.len() != f.tail.len() + length + rt.len()
                || !sequence.starts_with(f.tail.as_bytes())
                || !sequence.ends_with(rt.as_bytes())
                || !(0..length)
                    .all(|i| sequence[f.tail.len() + i] == source_bases[(g.forward_start + i) % n])
            {
                continue;
            }
            alternatives += 1;
            push(
                &mut overlay,
                n,
                g.forward_start,
                length,
                false,
                Role::Amplified,
                &format!("Selected exact product · role assignment {alternatives}"),
            );
            push(
                &mut overlay,
                n,
                g.forward_start,
                fl,
                false,
                Role::Forward,
                &label(fi, f),
            );
            push(
                &mut overlay,
                n,
                g.reverse_start,
                rl,
                true,
                Role::Reverse,
                &label(ri, r),
            );
            tail(&mut overlay, f, fi, "forward-role");
            tail(&mut overlay, r, ri, "reverse-role");
        }
    }
    if alternatives == 0 {
        return Err("Selected exact product has no compatible supplied-oligo role assignment at these coordinates.".into());
    }
    overlay.summary = format!(
        "Selected exact-PCR reference product · {alternatives} compatible directed oligo-role assignment(s). Other products remain separate."
    );
    overlay.warnings.push("Placement comes from the selected native result, not a new scan or evidence of experimental amplification/specificity. 5′ tails are not template bases.".into());
    if alternatives > 1 {
        overlay.warnings.push("The legacy exact result does not retain oligo-role indices. All compatible role assignments are shown; names are not unique identifiers and no assignment is silently preferred.".into());
    }
    finish(overlay, n)
}

fn mapped_location_matches(hit: &BindingHit, n: usize) -> bool {
    let location = if hit.reverse {
        let Location::Complement(inner) = &hit.location else {
            return false;
        };
        inner.as_ref()
    } else {
        &hit.location
    };
    let first = hit.length.min(n - hit.start);
    let range_matches = |loc: &Location, start, end| matches!(loc,Location::Range{start:a,end:b,before:false,after:false} if *a==start && *b==end);
    if first == hit.length {
        range_matches(location, hit.start, hit.start + first)
    } else {
        match location {
            Location::Join(parts) if parts.len() == 2 => {
                range_matches(&parts[0], hit.start, n)
                    && range_matches(&parts[1], 0, hit.length - first)
            }
            _ => false,
        }
    }
}
fn compare_domain(
    m: &Molecule,
    p: &Primer,
    start: usize,
    reverse: bool,
) -> Result<(Vec<usize>, Vec<usize>, usize), String> {
    domain_guard(m, start, p.annealing.len())?;
    let n = m.sequence.len();
    let length = p.annealing.len();
    let mut mismatches = Vec::new();
    let mut ambiguous = Vec::new();
    let mut exact = vec![false; length];
    for (i, base) in p.annealing.as_str().bytes().enumerate() {
        let pos = (start + if reverse { length - 1 - i } else { i }) % n;
        let reference = m.sequence.as_str().as_bytes()[pos];
        let reference = if reverse {
            collagene_model::complement(reference, Alphabet::Dna)
        } else {
            reference
        };
        let (a, b) = (base_mask(reference), base_mask(base));
        if a & b == 0 {
            mismatches.push(i);
        } else if a.count_ones() != 1 || b.count_ones() != 1 {
            ambiguous.push(i);
        } else {
            exact[i] = true;
        }
    }
    let three = exact.iter().rev().take_while(|&&v| v).count();
    Ok((mismatches, ambiguous, three))
}
fn markers(
    overlay: &mut Overlay,
    n: usize,
    start: usize,
    length: usize,
    reverse: bool,
    positions: &[usize],
    role: Role,
    label: &str,
) {
    for &i in positions {
        let pos = (start + if reverse { length - 1 - i } else { i }) % n;
        push(overlay, n, pos, 1, reverse, role, label);
    }
}
fn binding(
    overlay: &mut Overlay,
    m: &Molecule,
    p: &Primer,
    index: usize,
    hit: &BindingHit,
    direction: Role,
) -> Result<(), String> {
    domain_guard(m, hit.start, hit.length)?;
    if hit.length != p.annealing.len()
        || hit.mismatches.len() > 200
        || hit.ambiguous_positions.len() > 200
        || hit.reverse != (direction == Role::Reverse)
        || hit.crosses_origin != (hit.start + hit.length > m.sequence.len())
        || !mapped_location_matches(hit, m.sequence.len())
    {
        return Err("Selected binding geometry is inconsistent.".into());
    }
    let (mismatches, ambiguous, three) = compare_domain(m, p, hit.start, hit.reverse)?;
    if mismatches != hit.mismatches
        || ambiguous != hit.ambiguous_positions
        || three != hit.three_prime_exact
        || hit.certain != (mismatches.is_empty() && ambiguous.is_empty())
    {
        return Err("Selected binding evidence differs from actual source/oligo calls.".into());
    }
    let text = label(index, p);
    push(
        overlay,
        m.sequence.len(),
        hit.start,
        hit.length,
        hit.reverse,
        direction,
        &text,
    );
    markers(
        overlay,
        m.sequence.len(),
        hit.start,
        hit.length,
        hit.reverse,
        &mismatches,
        Role::Mismatch,
        &format!("{text} · disjoint primer/reference calls"),
    );
    markers(
        overlay,
        m.sequence.len(),
        hit.start,
        hit.length,
        hit.reverse,
        &ambiguous,
        Role::Ambiguous,
        &format!("{text} · compatible ambiguous calls"),
    );
    Ok(())
}
pub fn mismatch(
    m: &Molecule,
    source: &EventBinding,
    epoch: u64,
    forward: &Primer,
    reverse: &Primer,
    report: &PcrReport,
    index: usize,
) -> Result<Overlay, String> {
    primer_guard(forward)?;
    primer_guard(reverse)?;
    if report.source_revision != source.revision
        || report.products.len() > 10_000
        || report.warnings.len() > 32
        || report.warnings.iter().any(|s| s.len() > 4096)
        || report.bindings.iter().any(|b| {
            b.hits.len() > 10_000 || b.annealing.len() > 200 || b.source_revision != source.revision
        })
    {
        return Err("Mismatch-PCR report exceeds bounds or belongs to another source.".into());
    }
    let g = report
        .products
        .get(index)
        .ok_or("Select a reported mismatch-PCR product.")?;
    if g.forward_primer_index > 1
        || g.reverse_primer_index > 1
        || g.substitutions.len() > 400
        || [&g.forward_binding, &g.reverse_binding].iter().any(|hit| {
            hit.length > 200 || hit.mismatches.len() > 200 || hit.ambiguous_positions.len() > 200
        })
    {
        return Err("Selected oligo indices/change count are invalid.".into());
    }
    source_guard(m, source)?;
    let primers = [forward, reverse];
    let f = primers[g.forward_primer_index];
    let r = primers[g.reverse_primer_index];
    for (binding, p) in report.bindings.iter().zip(primers) {
        if binding.annealing != p.annealing.as_str() {
            return Err("Mismatch-PCR primer pair changed.".into());
        }
    }
    let n = m.sequence.len();
    domain_guard(m, g.forward_start, g.template_length)?;
    domain_guard(m, g.forward_start, f.annealing.len())?;
    domain_guard(m, g.reverse_start, r.annealing.len())?;
    if g.forward_binding.start != g.forward_start
        || g.reverse_binding.start != g.reverse_start
        || g.forward_tail_length != f.tail.len()
        || g.reverse_tail_length != r.tail.len()
        || g.sequence.alphabet() != Alphabet::Dna
        || g.sequence.len() != f.tail.len() + g.template_length + r.tail.len()
        || g.crosses_origin != (g.forward_start + g.template_length > n)
        || template_length(m, g.forward_start, g.reverse_start + r.annealing.len())
            != Some(g.template_length)
        || g.template_length < f.annealing.len().max(r.annealing.len())
    {
        return Err("Selected mismatch-PCR product geometry differs.".into());
    }
    let rd = rc(r.annealing.as_str())?;
    let rt = rc(&r.tail)?;
    let bases = g.sequence.as_str().as_bytes();
    if !bases.starts_with(f.tail.as_bytes()) || !bases.ends_with(rt.as_bytes()) {
        return Err("Selected product tails differ.".into());
    }
    let right = g.template_length - r.annealing.len();
    for j in 0..g.template_length {
        let forward_call = (j < f.annealing.len()).then(|| f.annealing.as_str().as_bytes()[j]);
        let reverse_call = (j >= right).then(|| rd.as_bytes()[j - right]);
        if forward_call.zip(reverse_call).is_some_and(|(a, b)| a != b) {
            return Err("Conflicting primer overlap cannot be projected.".into());
        }
        let expected = forward_call
            .or(reverse_call)
            .unwrap_or(m.sequence.as_str().as_bytes()[(g.forward_start + j) % n]);
        if bases[f.tail.len() + j] != expected {
            return Err("Product calls differ from the selected source/primer geometry.".into());
        }
    }
    let mut overlay = new(source, epoch, Kind::MismatchProduct);
    push(
        &mut overlay,
        n,
        g.forward_start,
        g.template_length,
        false,
        Role::Amplified,
        "Selected mismatch-PCR template footprint",
    );
    binding(
        &mut overlay,
        m,
        f,
        g.forward_primer_index,
        &g.forward_binding,
        Role::Forward,
    )?;
    binding(
        &mut overlay,
        m,
        r,
        g.reverse_primer_index,
        &g.reverse_binding,
        Role::Reverse,
    )?;
    tail(&mut overlay, f, g.forward_primer_index, "forward-role");
    tail(&mut overlay, r, g.reverse_primer_index, "reverse-role");
    overlay.summary = format!(
        "Selected mismatch-PCR product {} · search complete: {} · concrete enumeration complete: {}.",
        index + 1,
        report.search_complete,
        report.enumeration_complete
    );
    overlay.warnings = report.warnings.clone();
    overlay.warnings.push("Markers use each annealing domain's own 5′→3′ indices mapped to original source coordinates. Compatibility is not extension probability or experimental amplification.".into());
    if !report.search_complete || !report.enumeration_complete {
        overlay.warnings.push("This selected geometry cannot establish unique amplification or absence of other products.".into());
    }
    finish(overlay, n)
}

pub fn captured(
    m: &Molecule,
    source: &EventBinding,
    epoch: u64,
    draft: &PairDraft,
    forward: &Primer,
    reverse: &Primer,
) -> Result<Option<Overlay>, String> {
    primer_guard(forward)?;
    primer_guard(reverse)?;
    let supplied = [forward, reverse];
    let drafts = [&draft.forward, &draft.reverse];
    for (d, p) in drafts.iter().zip(supplied) {
        if let Some(anchor) = &d.anchor {
            crate::primers::primer_pair_model::selected_guard(
                &anchor.source,
                anchor.start,
                anchor.end,
                &anchor.reference,
            )?;
        }
        let applied = d.primer()?;
        if applied.name != p.name || applied.annealing != p.annealing || applied.tail != p.tail {
            return Ok(None);
        }
    }
    source_guard(m, source)?;
    let mut overlay = new(source, epoch, Kind::CapturedSelection);
    for (i, (d, p)) in drafts.iter().zip(supplied).enumerate() {
        let Some(anchor) = &d.anchor else {
            continue;
        };
        crate::primers::primer_pair_model::selected_guard(
            &anchor.source,
            anchor.start,
            anchor.end,
            &anchor.reference,
        )?;
        if &anchor.source != source {
            return Err("Captured selection belongs to a different source event.".into());
        }
        let oriented = crate::primers::primer_pair_model::anchor_reference(
            m,
            source,
            anchor.start,
            anchor.end,
            anchor.reverse,
        )?;
        if oriented != anchor.reference || p.annealing.len() != anchor.reference.len() {
            return Err("Captured source calls or annealing length changed; select the intended domain again.".into());
        }
        let (different, ambiguous, _) = compare_domain(m, p, anchor.start, anchor.reverse)?;
        let text = format!("{} · captured selection comparison only", label(i, p));
        push(
            &mut overlay,
            m.sequence.len(),
            anchor.start,
            p.annealing.len(),
            anchor.reverse,
            if anchor.reverse {
                Role::Reverse
            } else {
                Role::Forward
            },
            &text,
        );
        markers(
            &mut overlay,
            m.sequence.len(),
            anchor.start,
            p.annealing.len(),
            anchor.reverse,
            &different,
            Role::Mismatch,
            &text,
        );
        markers(
            &mut overlay,
            m.sequence.len(),
            anchor.start,
            p.annealing.len(),
            anchor.reverse,
            &ambiguous,
            Role::Ambiguous,
            &text,
        );
        tail(&mut overlay, p, i, "captured oligo");
    }
    if overlay.spans.is_empty() {
        return Ok(None);
    }
    overlay.summary =
        "Captured source-selection comparison · no binding search or amplification assessment."
            .into();
    overlay.warnings.push("These anchors describe the selected reference regions only. Difference/ambiguity markers compare current applied oligos with those regions; they do not establish primer binding, extension, specificity or an amplified interval.".into());
    finish(overlay, m.sequence.len()).map(Some)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::primers::primer_pair_model::SelectionSeed;
    use collagene_core::{
        cloning::pcr::{mutagenesis, primers::BindingOptions},
        documents::workspace::Workspace,
        verification,
    };
    fn fixture(bases: &str, topology: Topology) -> (Molecule, EventBinding) {
        let m = Molecule::dna("overlay source", bases, topology).unwrap();
        let w = Workspace::from_documents("overlay", vec![m.clone()], 0).unwrap();
        let d = &w.documents[0];
        let binding = EventBinding {
            document_id: d.id.clone(),
            event_id: d.head_event.clone().unwrap(),
            revision: m.revision(),
        };
        (m, binding)
    }
    fn p(name: &str, domain: &str, tail: &str) -> Primer {
        Primer::new(name, domain, tail).unwrap()
    }
    fn mismatch_fixture(bases: &str) -> (Molecule, EventBinding, Primer, Primer, PcrReport, usize) {
        let (m, b) = fixture(bases, Topology::Linear);
        let (f, r) = (p("alias", "ATGA", "CGT"), p("alias", "GTCA", "GAC"));
        let report = mutagenesis::simulate_pcr(
            &m,
            &f,
            &r,
            &mutagenesis::PcrOptions {
                binding: BindingOptions {
                    max_mismatches: 1,
                    min_three_prime_exact: 1,
                    ..Default::default()
                },
                ..Default::default()
            },
        )
        .unwrap();
        let index = report
            .products
            .iter()
            .position(|g| {
                g.forward_start == 0
                    && g.reverse_start == 12
                    && g.forward_primer_index == 0
                    && g.reverse_primer_index == 1
            })
            .unwrap();
        (m, b, f, r, report, index)
    }
    #[test]
    fn exact_selected_geometry_maps_both_domains_and_excludes_tails_from_source() {
        let (m, b) = fixture("ACGACCCCGGGGTAAC", Topology::Linear);
        let (f, r) = (p("same name", "ACGA", "CGT"), p("same name", "GTTA", "GAC"));
        let products = verification::pcr::simulate_pcr_documents(
            &m,
            &f,
            &r,
            collagene_core::chemistry::Strand::Top,
        )
        .unwrap();
        let selected = products
            .iter()
            .find(|g| g.geometry.forward_start == 0 && g.geometry.reverse_start == 12)
            .unwrap();
        let overlay = exact(&m, &b, 7, &f, &r, selected).unwrap();
        assert_eq!(overlay.kind, Kind::ExactProduct);
        assert!(
            overlay
                .spans
                .iter()
                .any(|s| s.role == Role::Forward && s.start == 0 && s.end == 4 && !s.reverse)
        );
        assert!(
            overlay
                .spans
                .iter()
                .any(|s| s.role == Role::Reverse && s.start == 12 && s.end == 16 && s.reverse)
        );
        assert!(
            overlay
                .spans
                .iter()
                .any(|s| s.role == Role::Amplified && s.start == 0 && s.end == 16)
        );
        assert_eq!(overlay.tails.len(), 2);
        assert!(overlay.tails.iter().any(|tail| tail.ends_with("CGT")));
        assert!(overlay.tails.iter().any(|tail| tail.ends_with("GAC")));
        assert!(overlay.validate_for_length(16));
        assert!(overlay.is_current(&b, 7));
        assert!(!overlay.is_current(&b, 8));
    }
    #[test]
    fn same_oligos_keep_all_four_compatible_roles_and_origin_spans() {
        let (m, b) = fixture("CCCTCGTACGA", Topology::Circular);
        let (f, r) = (p("first", "ACGA", "GG"), p("second", "ACGA", "GG"));
        let products = verification::pcr::simulate_pcr_documents(
            &m,
            &f,
            &r,
            collagene_core::chemistry::Strand::Top,
        )
        .unwrap();
        let selected = products
            .iter()
            .find(|g| g.geometry.forward_start == 7 && g.geometry.reverse_start == 3)
            .unwrap();
        assert!(selected.geometry.crosses_origin);
        let overlay = exact(&m, &b, 0, &f, &r, selected).unwrap();
        assert!(overlay.summary.contains("4 compatible"));
        assert!(
            overlay
                .warnings
                .iter()
                .any(|s| s.contains("All compatible role assignments"))
        );
        assert_eq!(
            overlay
                .spans
                .iter()
                .filter(|s| s.role == Role::Forward)
                .count(),
            4
        );
        assert!(
            overlay
                .spans
                .iter()
                .any(|s| s.role == Role::Amplified && s.start == 7 && s.end == 11)
        );
        assert!(
            overlay
                .spans
                .iter()
                .any(|s| s.role == Role::Amplified && s.start == 0 && s.end == 7)
        );
    }
    #[test]
    fn mismatch_indices_map_reverse_oligo_five_prime_coordinates_to_source() {
        let (m, b, f, r, report, index) = mismatch_fixture("ACGACCCCGGGGTAAC");
        assert_eq!(report.products[index].reverse_binding.mismatches, vec![2]);
        let overlay = mismatch(&m, &b, 9, &f, &r, &report, index).unwrap();
        let markers = overlay
            .spans
            .iter()
            .filter(|s| s.role == Role::Mismatch)
            .map(|s| (s.start, s.end, s.reverse))
            .collect::<Vec<_>>();
        assert_eq!(markers, vec![(1, 2, false), (13, 14, true)]);
        assert!(
            overlay
                .spans
                .iter()
                .any(|s| s.role == Role::Forward && s.label.starts_with("Oligo 1"))
        );
        assert!(
            overlay
                .spans
                .iter()
                .any(|s| s.role == Role::Reverse && s.label.starts_with("Oligo 2"))
        );
        assert_eq!(overlay.kind, Kind::MismatchProduct);
    }
    #[test]
    fn ambiguous_calls_and_partial_enumeration_remain_explicit() {
        let (m, b, f, r, mut report, index) = mismatch_fixture("ANGACCCCGGGGTAAC");
        report.search_complete = false;
        report.enumeration_complete = false;
        let overlay = mismatch(&m, &b, 0, &f, &r, &report, index).unwrap();
        assert!(
            overlay
                .spans
                .iter()
                .any(|s| s.role == Role::Ambiguous && s.start == 1 && s.end == 2)
        );
        assert!(overlay.summary.contains("search complete: false"));
        assert!(overlay.summary.contains("enumeration complete: false"));
        assert!(
            overlay
                .warnings
                .iter()
                .any(|s| s.contains("cannot establish unique"))
        );
    }
    #[test]
    fn captured_anchors_are_comparisons_never_an_amplified_or_binding_result() {
        let (m, b) = fixture("ACGACCCCGGGGTAAC", Topology::Linear);
        let mut draft = PairDraft {
            forward: SelectionSeed::capture(&m, &b, 0, 4)
                .unwrap()
                .oligo("forward", false, "CGT")
                .unwrap(),
            reverse: SelectionSeed::capture(&m, &b, 12, 16)
                .unwrap()
                .oligo("reverse", true, "GAC")
                .unwrap(),
        };
        draft.forward.annealing = "ATGA".into();
        draft.reverse.annealing = "GTCA".into();
        let (f, r) = draft.primers().unwrap();
        let overlay = captured(&m, &b, 5, &draft, &f, &r).unwrap().unwrap();
        assert_eq!(overlay.kind, Kind::CapturedSelection);
        assert!(!overlay.spans.iter().any(|s| s.role == Role::Amplified));
        assert!(overlay.summary.contains("no binding search"));
        assert_eq!(
            overlay
                .spans
                .iter()
                .filter(|s| s.role == Role::Mismatch)
                .map(|s| s.start)
                .collect::<Vec<_>>(),
            vec![1, 13]
        );
        let changed = p("forward", "ACGA", "CGT");
        assert!(captured(&m, &b, 5, &draft, &changed, &r).unwrap().is_none());
        let mut stale = draft.clone();
        stale.reverse.anchor.as_mut().unwrap().source.event_id =
            format!("event-{}", "0".repeat(64));
        assert!(captured(&m, &b, 5, &stale, &f, &r).is_err());
        stale = draft.clone();
        stale.reverse.anchor.as_mut().unwrap().reference = "AAAA".into();
        assert!(captured(&m, &b, 5, &stale, &f, &r).is_err());
    }
    #[test]
    fn changed_source_and_invalid_selected_evidence_cannot_publish() {
        let (m, b, f, r, mut report, index) = mismatch_fixture("ACGACCCCGGGGTAAC");
        report.products[index].reverse_binding.mismatches = vec![0];
        assert!(mismatch(&m, &b, 0, &f, &r, &report, index).is_err());
        let (m, b, f, r, mut report, index) = mismatch_fixture("ACGACCCCGGGGTAAC");
        report.products[index].reverse_start = usize::MAX;
        assert!(mismatch(&m, &b, 0, &f, &r, &report, index).is_err());
        let mut changed = b.clone();
        changed.revision = "0".repeat(64);
        assert!(mismatch(&m, &changed, 0, &f, &r, &report, index).is_err());
        assert!(mismatch(&m, &b, 0, &f, &r, &report, usize::MAX).is_err());
    }
    #[test]
    fn context_guard_uses_event_identity_and_render_validation_is_bounded() {
        let (m, b) = fixture("ACGACCCCGGGGTAAC", Topology::Linear);
        let mut overlay = new(&b, 1, Kind::CapturedSelection);
        push(&mut overlay, 16, 0, 4, false, Role::Forward, "captured");
        let mut alias = b.clone();
        alias.document_id = "other_document".into();
        assert!(!overlay.is_current(&alias, 1));
        assert!(overlay.validate_for_length(m.sequence.len()));
        overlay.spans[0].role = Role::Amplified;
        assert!(!overlay.validate_for_length(16));
        overlay.spans[0].role = Role::Forward;
        overlay.spans[0].end = 17;
        assert!(!overlay.validate_for_length(16));
        overlay.spans = vec![
            Span {
                start: 0,
                end: 1,
                reverse: false,
                role: Role::Forward,
                label: String::new()
            };
            MAX_SPANS + 1
        ];
        assert!(!overlay.validate_for_length(16));
    }

    #[test]
    fn captured_wrapped_domains_split_at_origin_and_map_reverse_mismatches_modularly() {
        let (m, b) = fixture("ACGTAAGC", Topology::Circular);
        let seed = SelectionSeed::capture(&m, &b, 6, 11).unwrap();
        let mut draft = PairDraft {
            forward: seed.oligo("forward", false, "TT").unwrap(),
            reverse: seed.oligo("reverse", true, "GG").unwrap(),
        };
        draft.forward.annealing = "GTACG".into();
        draft.reverse.annealing = "CATGC".into();
        let (f, r) = draft.primers().unwrap();
        let overlay = captured(&m, &b, 2, &draft, &f, &r).unwrap().unwrap();
        assert_eq!(overlay.kind, Kind::CapturedSelection);
        assert!(
            overlay
                .spans
                .iter()
                .all(|span| span.role != Role::Amplified)
        );
        assert_eq!(
            overlay
                .spans
                .iter()
                .filter(|span| span.role == Role::Forward)
                .map(|span| (span.start, span.end, span.reverse))
                .collect::<Vec<_>>(),
            vec![(6, 8, false), (0, 3, false)]
        );
        assert_eq!(
            overlay
                .spans
                .iter()
                .filter(|span| span.role == Role::Reverse)
                .map(|span| (span.start, span.end, span.reverse))
                .collect::<Vec<_>>(),
            vec![(6, 8, true), (0, 3, true)]
        );
        assert_eq!(
            overlay
                .spans
                .iter()
                .filter(|span| span.role == Role::Mismatch)
                .map(|span| (span.start, span.reverse))
                .collect::<Vec<_>>(),
            vec![(7, false), (1, true)]
        );
        assert!(overlay.validate_for_length(8));
    }
}
