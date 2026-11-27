//! Explicit ideal ligase sealing of true duplex nicks after fresh TOPO transfer.
//!
//! Chemical prerequisite: juxtaposed3′OH and5′phosphate in duplex DNA. See NEB's
//! primary ligase reference, https://www.neb.com/en/tools-and-resources/selection-charts/t4-dna-ligase-products-table
//! and Lohman et al., https://pmc.ncbi.nlm.nih.gov/articles/PMC3243518/ .
//! This narrower exact-template profile does not model ATP turnover, efficiency,
//! mismatched/gapped substrates, repair, phosphorylation or kit performance.
use crate::{
    Budget, Call, ChainEnds, Error, FivePrime, Graph, Limits, Outcome, ProductGraph, Substrates,
    ThreePrime, Transfer,
};
use serde::{Deserialize, Serialize};
use std::collections::BTreeSet;
pub const ENGINE: &str = "topo-explicit-nick-sealing-1";
pub const MODEL: &str = "paired-continuous-template-three-oh-five-p-1";
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,PartialEq,Eq,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct NickSelection{pub three:Call,pub five:Call}}
dto! {pub struct Request{pub ligase_assumption:String,pub nicks:Vec<NickSelection>}}
dto! {pub struct SealEvidence{pub step:usize,pub selection:NickSelection,pub outcome:Outcome,pub reason:String,
pub three_prime:ThreePrime,pub five_prime:FivePrime,pub template_three:Option<Call>,pub template_five:Option<Call>}}
dto! {pub struct Report{pub engine:String,pub model:String,pub context_revision:String,pub parent:crate::Report,pub request:Request,
pub steps:Vec<SealEvidence>,pub product:ProductGraph,pub source_bases:usize,
pub all_requested_nicks_sealed:bool,pub unsupported_or_unknown:bool,pub assumptions:Vec<String>}}
fn bad(s: &str) -> Error {
    crate::bad(&format!("nick sealing: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn valid_call(call: Call, l: &Limits) -> bool {
    call.strand < l.max_strands && call.index < l.max_source_bases
}
/// Cheap borrowed request guard before source hashes/preparation. Repeated nick
/// requests are retained in order; a consumed endpoint blocks the later step.
pub fn request_preflight(r: &Request, l: &Limits) -> Result<(), Error> {
    crate::limits(l)?;
    if !crate::text(&r.ligase_assumption, 4096) || r.ligase_assumption.chars().any(char::is_control)
    {
        return Err(bad(
            "provide an explicit bounded ligase/conditions assumption",
        ));
    }
    crate::bounded(r.nicks.len(), l.max_steps, Budget::Steps)?;
    if r.nicks
        .iter()
        .any(|n| !valid_call(n.three, l) || !valid_call(n.five, l))
    {
        return Err(bad("nick call exceeds declared source bounds"));
    }
    Ok(())
}
/// Combined reservations cover parent simulation, retained parent output, fresh
/// graph reconstruction, sealing, child output and both ancestry audits.
pub fn preflight(
    s: &Substrates,
    original: &[Transfer],
    l: &Limits,
    r: &Request,
) -> Result<usize, Error> {
    request_preflight(r, l)?;
    let n = crate::preflight(s, original, l)?;
    let steps = crate::add(original.len(), r.nicks.len())?;
    crate::bounded(steps, l.max_steps, Budget::Steps)?;
    for selection in &r.nicks {
        for call in [selection.three, selection.five] {
            if s.strands
                .get(call.strand)
                .is_none_or(|strand| call.index >= strand.sequence.len())
            {
                return Err(bad(
                    "nick call is not an actual immutable source nucleotide",
                ));
            }
        }
    }
    crate::bounded(
        crate::add(crate::mul(n, 128)?, crate::mul(steps, 128)?)?,
        l.max_work_cells,
        Budget::WorkCells,
    )?;
    let bytes = crate::add(
        crate::add(crate::mul(n, 256)?, crate::mul(s.pairing.len(), 256)?)?,
        crate::add(
            crate::add(crate::mul(s.strands.len(), 8192)?, crate::mul(steps, 4096)?)?,
            16384,
        )?,
    )?;
    crate::bounded(bytes, l.max_work_bytes, Budget::WorkBytes)?;
    crate::bounded(crate::mul(n, 2)?, l.max_output_bases, Budget::OutputBases)?;
    Ok(n)
}
fn context(parent: &str, l: &Limits, r: &Request) -> Result<String, Error> {
    crate::hash(&(ENGINE, MODEL, parent, l, r))
}
pub fn revision(
    s: &Substrates,
    original: &[Transfer],
    l: &Limits,
    r: &Request,
) -> Result<String, Error> {
    preflight(s, original, l, r)?;
    context(&crate::batch_revision(s, original, l)?, l, r)
}
fn chemistry(s: &Substrates, g: &Graph, selection: &NickSelection) -> (ThreePrime, FivePrime) {
    let three = crate::terminal3(s, g, g.flat(selection.three));
    let five = if selection.five.index == 0 {
        s.strands[selection.five.strand].five_prime
    } else {
        FivePrime::Unknown
    };
    (three, five)
}
fn eligible(s: &Substrates, g: &Graph, n: &NickSelection) -> Result<(), (Outcome, &'static str)> {
    let (tc, fc) = chemistry(s, g, n);
    eligible_with_chemistry(s, g, n, tc, fc)
}
pub(crate) fn eligible_with_chemistry(
    s: &Substrates,
    g: &Graph,
    n: &NickSelection,
    tc: ThreePrime,
    fc: FivePrime,
) -> Result<(), (Outcome, &'static str)> {
    let three = g.flat(n.three);
    let five = g.flat(n.five);
    if three == five {
        return Err((
            Outcome::Unsupported,
            "one nucleotide cannot provide both sides of this nick profile",
        ));
    }
    if g.next[three].is_some() || g.prev[five].is_some() {
        return Err((
            Outcome::BlockedUnderDeclaredModel,
            "requested termini are not exposed or have already formed a covalent bond",
        ));
    }
    match tc {
        ThreePrime::Hydroxyl => {}
        ThreePrime::Unknown => {
            return Err((Outcome::Unknown, "actual exposed3′chemistry is unknown"));
        }
        ThreePrime::Phosphotyrosyl { .. } => {
            return Err((
                Outcome::Unsupported,
                "ligase sealing does not consume a TOPO activation; an explicit strand transfer is separate",
            ));
        }
        _ => {
            return Err((
                Outcome::BlockedUnderDeclaredModel,
                "nick sealing requires actual3′OH; no repair/dephosphorylation is inferred",
            ));
        }
    }
    match fc {
        FivePrime::Phosphate => {}
        FivePrime::Unknown => {
            return Err((Outcome::Unknown, "actual exposed5′chemistry is unknown"));
        }
        _ => {
            return Err((
                Outcome::BlockedUnderDeclaredModel,
                "nick sealing requires actual5′phosphate; no phosphorylation or unblocking is inferred",
            ));
        }
    }
    let Some(tp) = g.pairs[three] else {
        return Err((
            Outcome::Unsupported,
            "3′terminus is unpaired, not a duplex-templated nick",
        ));
    };
    let Some(fp) = g.pairs[five] else {
        return Err((
            Outcome::Unsupported,
            "5′terminus is unpaired, not a duplex-templated nick",
        ));
    };
    if tp == fp
        || g.pairs[tp] != Some(three)
        || g.pairs[fp] != Some(five)
        || g.next[fp] != Some(tp)
        || g.prev[tp] != Some(fp)
    {
        return Err((
            Outcome::Unsupported,
            "opposite strand is not covalently continuous in the antiparallel direction across this exact nick",
        ));
    }
    for (a, b) in [(three, tp), (five, fp)] {
        let Some(complement) = crate::comp(g.base(s, a)) else {
            return Err((
                Outcome::Unknown,
                "ambiguous nick/template calls do not establish exact duplex ligation geometry",
            ));
        };
        if crate::comp(g.base(s, b)).is_none() {
            return Err((Outcome::Unknown, "ambiguous opposite-template call"));
        }
        if complement != g.base(s, b) {
            return Err((
                Outcome::Unsupported,
                "nick/template calls are not exact Watson–Crick pairs",
            ));
        }
    }
    Ok(())
}
/// No serialized graph/report input exists. Parent evidence is generated in this
/// call, then only its actually applied changes rebuild the private source graph.
pub fn simulate(
    s: &Substrates,
    original: &[Transfer],
    l: &Limits,
    r: &Request,
) -> Result<Report, Error> {
    let fresh = crate::simulate_with_checkpoints(s, original, l, &[])?;
    seal_fresh(s, original, l, r, fresh)
}
pub fn seal_fresh(
    s: &Substrates,
    original: &[Transfer],
    l: &Limits,
    r: &Request,
    mut fresh: crate::Executed,
) -> Result<Report, Error> {
    let n = preflight(s, original, l, r)?;
    if fresh.report().batch_revision != crate::batch_revision(s, original, l)? {
        return Err(bad("fresh native nick sealing parent differs"));
    }
    fresh.reserve_work(
        crate::add(
            crate::mul(n, 256)?,
            crate::mul(original.len() + r.nicks.len(), 256)?,
        )?,
        crate::mul(n, 512)?,
    )?;
    let parent = fresh.into_report();
    let context_revision = context(&parent.batch_revision, l, r)?;
    let mut g = Graph::new(s, n)?;
    for step in parent
        .steps
        .iter()
        .filter(|step| step.outcome == Outcome::Applied)
    {
        for pair in &step.removed_pairs {
            let a = g.flat(pair.a);
            let b = g.flat(pair.b);
            g.pairs[a] = None;
            g.pairs[b] = None;
            g.displaced[b] = Some(step.step);
        }
        for pair in &step.added_pairs {
            let a = g.flat(pair.a);
            let b = g.flat(pair.b);
            g.pairs[a] = Some(b);
            g.pairs[b] = Some(a);
        }
        let d = g.flat(step.donor_call);
        let a = g.flat(step.acceptor_call);
        g.next[d] = Some(a);
        g.prev[a] = Some(d);
    }
    let mut steps = Vec::with_capacity(r.nicks.len());
    for (step, selection) in r.nicks.iter().enumerate() {
        let (three_prime, five_prime) = chemistry(s, &g, selection);
        let d = g.flat(selection.three);
        let a = g.flat(selection.five);
        let template_three = g.pairs[d].map(|p| g.call(p));
        let template_five = g.pairs[a].map(|p| g.call(p));
        let (outcome, reason) = match eligible(s, &g, selection) {
            Ok(()) => {
                g.next[d] = Some(a);
                g.prev[a] = Some(d);
                (
                    Outcome::Applied,
                    "one declared3′OH→5′phosphate nick sealed; pairing, calls, activation and displaced material unchanged",
                )
            }
            Err(result) => result,
        };
        steps.push(SealEvidence {
            step,
            selection: selection.clone(),
            outcome,
            reason: reason.into(),
            three_prime,
            five_prime,
            template_three,
            template_five,
        });
    }
    let product = crate::output(s, &g, parent.product.released_activations.clone())?;
    let all_requested_nicks_sealed =
        !steps.is_empty() && steps.iter().all(|step| step.outcome == Outcome::Applied);
    let unsupported_or_unknown = parent.unsupported_or_unknown
        || steps
            .iter()
            .any(|step| matches!(step.outcome, Outcome::Unknown | Outcome::Unsupported));
    let report=Report{engine:ENGINE.into(),model:MODEL.into(),context_revision,parent,request:r.clone(),steps,product,source_bases:n,all_requested_nicks_sealed,unsupported_or_unknown,
        assumptions:vec!["The supplied ligase/conditions are an explicit operator hypothesis. This exact nick model assumes active ligase/cofactor availability; no ATP/AMP turnover, yield, kinetic mixture or experimental success is established.".into(),"Only exposed3′OH and5′phosphate calls paired to an actual continuous antiparallel DNA template are sealed. No phosphorylation, phosphotyrosyl activation consumption, repair, polymerase extension or mismatch/gap/flap removal is inferred.".into(),"All original source calls, parent TOPO transfer outcomes, pairing changes, displaced/unpaired material and released activation identities remain retained. Covalent edges alone change backbone topology.".into(),"Unknown, unsupported and blocked operations leave their incoming graph unchanged. Later selections describe this retained-state branch, not a distribution of possible experimental outcomes.".into()]};
    validate_report(s, original, l, r, &report)?;
    Ok(report)
}
fn three_shape(value: &ThreePrime) -> Result<(), Error> {
    if let ThreePrime::Phosphotyrosyl { activation_id, .. } = value {
        if !crate::text(activation_id, 256) {
            return Err(bad("activation text bound"));
        }
    }
    Ok(())
}
fn product_shape(product: &ProductGraph, l: &Limits, n: usize) -> Result<(), Error> {
    crate::bounded(product.chains.len(), l.max_strands, Budget::Strands)?;
    crate::bounded(
        product.source_pairing.len(),
        l.max_pair_spans,
        Budget::PairSpans,
    )?;
    crate::bounded(product.nicks.len(), l.max_strands, Budget::Strands)?;
    crate::bounded(product.unpaired.len(), n, Budget::SourceBases)?;
    crate::bounded(product.components.len(), l.max_strands, Budget::Strands)?;
    let mut total = 0;
    for chain in &product.chains {
        total = crate::add(total, chain.sequence.len())?;
        crate::bounded(total, n, Budget::SourceBases)?;
        if chain.sequence.alphabet() != collagene_model::Alphabet::Dna
            || chain.sources.len() > l.max_strands
        {
            return Err(bad("chain DNA/source placement shape"));
        }
        if let ChainEnds::Linear { three_prime, .. } = &chain.ends {
            three_shape(three_prime)?;
        }
    }
    if total != n {
        return Err(bad("product source-call count changed"));
    }
    for span in &product.source_pairing {
        if span.strand_a >= l.max_strands
            || span.strand_b >= l.max_strands
            || span.length == 0
            || span.start_a.checked_add(span.length).is_none_or(|v| v > n)
            || span.start_b.checked_add(span.length).is_none_or(|v| v > n)
        {
            return Err(bad("pair span shape changed"));
        }
    }
    for nick in &product.nicks {
        if !valid_call(nick.three, l) || !valid_call(nick.five, l) {
            return Err(bad("nick call shape"));
        }
        three_shape(&nick.three_prime)?;
    }
    for component in &product.components {
        if component.chains.len() > product.chains.len()
            || component.chains.iter().any(|i| *i >= product.chains.len())
            || component
                .paired_nucleotides
                .checked_add(component.unpaired_nucleotides)
                .is_none_or(|v| v > n)
        {
            return Err(bad("component shape changed"));
        }
    }
    for run in &product.unpaired {
        if run.chain >= product.chains.len()
            || run.start >= run.end
            || run.end > product.chains[run.chain].sequence.len()
            || run.displaced_by.len() > l.max_steps
            || run.displaced_by.iter().any(|i| *i >= l.max_steps)
        {
            return Err(bad("unpaired/displaced call run shape changed"));
        }
    }
    if product.released_activations.len() > l.max_steps
        || product
            .released_activations
            .iter()
            .any(|s| !crate::text(s, 256))
    {
        return Err(bad("released activation shape"));
    }
    Ok(())
}
impl Report {
    /// Bounded readonly binding/ancestry-shape check, not chemical replay. Only
    /// simulate/replay reconstruct from actual original substrates and selections.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        let l = &self.parent.limits;
        request_preflight(&self.request, l)?;
        crate::bounded(self.source_bases, l.max_source_bases, Budget::SourceBases)?;
        crate::bounded(
            crate::add(self.parent.steps.len(), self.steps.len())?,
            l.max_steps,
            Budget::Steps,
        )?;
        if self.engine != ENGINE
            || self.model != MODEL
            || !hex(expected)
            || self.context_revision != expected
            || !hex(&self.parent.batch_revision)
            || self.context_revision != context(&self.parent.batch_revision, l, &self.request)?
            || self.parent.engine != crate::ENGINE
            || self.parent.model != crate::PRIMARY_MODEL
            || self.source_bases != self.parent.source_bases
            || self.steps.len() != self.request.nicks.len()
            || self.parent.source_revisions.len() > l.max_strands
            || self.parent.source_revisions.iter().any(|s| !hex(s))
            || self.parent.assumptions.len() > 16
            || self.parent.assumptions.iter().any(|s| s.len() > 4096)
            || self.assumptions.len() > 16
            || self.assumptions.iter().any(|s| s.len() > 4096)
        {
            return Err(bad("sealing source/model/context shape changed"));
        }
        for (index, step) in self.parent.steps.iter().enumerate() {
            if step.step != index
                || step.reason.len() > 4096
                || step.removed_pairs.len() > 64
                || step.added_pairs.len() > 64
                || !valid_call(step.donor_call, l)
                || !valid_call(step.acceptor_call, l)
            {
                return Err(bad("parent transfer evidence shape"));
            }
            for p in step.removed_pairs.iter().chain(&step.added_pairs) {
                if !valid_call(p.a, l) || !valid_call(p.b, l) {
                    return Err(bad("parent pair-change call shape"));
                }
            }
            if step
                .released_activation
                .as_ref()
                .is_some_and(|s| !crate::text(s, 256))
            {
                return Err(bad("activation ID shape"));
            }
            if let Some(d) = &step.displacement {
                if d.calls_five_to_three.len() > 64
                    || d.calls_five_to_three.iter().any(|c| !valid_call(*c, l))
                    || d.five_prime_continuation.is_some_and(|c| !valid_call(c, l))
                    || d.continuation_pair.is_some_and(|c| !valid_call(c, l))
                {
                    return Err(bad("parent displacement shape"));
                }
            }
        }
        for (i, step) in self.steps.iter().enumerate() {
            if step.step != i
                || step.selection != self.request.nicks[i]
                || step.reason.len() > 4096
                || step.template_three.is_some_and(|c| !valid_call(c, l))
                || step.template_five.is_some_and(|c| !valid_call(c, l))
            {
                return Err(bad("nick outcome identity/shape"));
            }
            three_shape(&step.three_prime)?;
            if step.outcome == Outcome::Applied
                && (step.three_prime != ThreePrime::Hydroxyl
                    || step.five_prime != FivePrime::Phosphate
                    || step.template_three.is_none()
                    || step.template_five.is_none())
            {
                return Err(bad("applied nick has unsupported chemical evidence"));
            }
        }
        if self.all_requested_nicks_sealed
            != (!self.steps.is_empty() && self.steps.iter().all(|s| s.outcome == Outcome::Applied))
            || self.unsupported_or_unknown
                != (self.parent.unsupported_or_unknown
                    || self
                        .steps
                        .iter()
                        .any(|s| matches!(s.outcome, Outcome::Unknown | Outcome::Unsupported)))
            || self.product.source_pairing != self.parent.product.source_pairing
            || self.product.released_activations != self.parent.product.released_activations
        {
            return Err(bad(
                "sealing altered parent pairing/activation or outcome flags",
            ));
        }
        product_shape(&self.parent.product, l, self.source_bases)?;
        product_shape(&self.product, l, self.source_bases)?;
        let applied: BTreeSet<_> = self
            .steps
            .iter()
            .filter(|s| s.outcome == Outcome::Applied)
            .map(|s| (s.selection.three, s.selection.five))
            .collect();
        let retained: Vec<_> = self
            .parent
            .product
            .nicks
            .iter()
            .filter(|n| !applied.contains(&(n.three, n.five)))
            .collect();
        if retained != self.product.nicks.iter().collect::<Vec<_>>()
            || applied.iter().any(|(three, five)| {
                !self
                    .parent
                    .product
                    .nicks
                    .iter()
                    .any(|n| n.three == *three && n.five == *five)
            })
        {
            return Err(bad(
                "sealing must remove exactly the applied original nicks",
            ));
        }
        crate::report_size(self, l.max_report_bytes)?;
        Ok(())
    }
}
fn audit_product(s: &Substrates, product: &ProductGraph, n: usize) -> Result<(), Error> {
    let mut coverage = vec![false; s.strands.len()];
    let mut total = 0;
    for chain in &product.chains {
        let mut cursor = 0;
        for p in &chain.sources {
            if p.source_strand >= s.strands.len()
                || coverage[p.source_strand]
                || p.start != 0
                || p.end != s.strands[p.source_strand].sequence.len()
                || p.output_start != cursor
            {
                return Err(bad(
                    "child source mapping is not exact once-only whole-strand ancestry",
                ));
            }
            let next = crate::add(cursor, p.end)?;
            if chain.sequence.as_str().as_bytes().get(cursor..next)
                != Some(s.strands[p.source_strand].sequence.as_str().as_bytes())
            {
                return Err(bad("child bases differ from original immutable calls"));
            }
            coverage[p.source_strand] = true;
            cursor = next;
        }
        if cursor != chain.sequence.len() {
            return Err(bad("child chain source mapping incomplete"));
        }
        total = crate::add(total, cursor)?;
    }
    if total != n || coverage.iter().any(|seen| !*seen) {
        return Err(bad("child source calls lost, duplicated or invented"));
    }
    Ok(())
}
/// Source-bound readonly guard, not fresh chemical eligibility. Used to install
/// only the response for the actual source/request token. Replay is separate.
pub fn validate_report(
    s: &Substrates,
    original: &[Transfer],
    l: &Limits,
    r: &Request,
    report: &Report,
) -> Result<(), Error> {
    let n = preflight(s, original, l, r)?;
    report.validate_binding(&revision(s, original, l, r)?)?;
    if &report.parent.limits != l
        || &report.request != r
        || report.parent.steps.len() != original.len()
        || report
            .parent
            .steps
            .iter()
            .zip(original)
            .any(|(a, b)| &a.transfer != b)
    {
        return Err(bad("parent source transfers/request changed"));
    }
    crate::audit_sources(s, &report.parent)?;
    audit_product(s, &report.product, n)?;
    Ok(())
}
pub fn replay(
    s: &Substrates,
    original: &[Transfer],
    l: &Limits,
    r: &Request,
    report: &Report,
) -> Result<(), Error> {
    validate_report(s, original, l, r, report)?;
    if &simulate(s, original, l, r)? != report {
        return Err(bad(
            "imported report differs from fresh native nick-sealing reconstruction",
        ));
    }
    Ok(())
}
