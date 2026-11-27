//! Ordered, source-conserving ideal processing after fresh elementary transfers.
//! Cleavage is an explicit hydrolytic hypothesis, not a named nuclease/kit claim.
use crate::{
    Budget, Call, ChainEnds, DisplacementKind, Error, FivePrime, Graph, Limits, Outcome,
    ProductGraph, Substrates, TerminalOverrides, ThreePrime, Transfer, nick_sealing as nick,
};
use serde::{Deserialize, Serialize};
pub const ENGINE: &str = "topo-source-processing-1";
pub const MODEL: &str = "ideal-anchored-three-flap-hydrolysis-and-templated-nick-sealing-1";
macro_rules! dto {($($t:tt)*)=>{#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Bond{pub upstream:Call,pub downstream:Call}}
dto! {#[derive(Copy)] pub enum FlapCleavageProfile{IdealDuplexAnchoredThreePrimeHydrolysisV1}}
dto! {pub enum Action{
    CleaveThreePrimeFlap{bond:Bond,displacement_step:usize,profile:FlapCleavageProfile,assumption:String},
    SealNick{selection:nick::NickSelection,ligase_assumption:String},
}}
dto! {pub struct ProcessingLimits{pub max_actions:usize,pub max_product_chains:usize,pub max_mapping_segments:usize,pub max_report_bytes:usize}}
impl Default for ProcessingLimits {
    fn default() -> Self {
        Self {
            max_actions: 16,
            max_product_chains: 64,
            max_mapping_segments: 4096,
            max_report_bytes: 32 * 1024 * 1024,
        }
    }
}
dto! {pub struct Request{pub actions:Vec<Action>,pub limits:ProcessingLimits}}
dto! {pub struct FlapProposal{pub displacement_step:usize,pub bond:Bond,pub flap_calls:Vec<Call>}}
dto! {pub struct FlapCleavageEvidence{
    pub step:usize,pub bond:Bond,pub displacement_step:usize,pub profile:FlapCleavageProfile,
    pub outcome:Outcome,pub reason:String,pub flap_calls:Vec<Call>,pub anchor_template:Option<Call>,
    pub created_three_prime:Option<ThreePrime>,pub created_five_prime:Option<FivePrime>,
}}
dto! {pub enum ActionEvidence{FlapCleavage(FlapCleavageEvidence),NickSealing(nick::SealEvidence)}}
impl ActionEvidence {
    pub fn outcome(&self) -> Outcome {
        match self {
            Self::FlapCleavage(e) => e.outcome,
            Self::NickSealing(e) => e.outcome,
        }
    }
    pub fn step(&self) -> usize {
        match self {
            Self::FlapCleavage(e) => e.step,
            Self::NickSealing(e) => e.step,
        }
    }
}
dto! {pub struct Report{
    pub engine:String,pub model:String,pub context_revision:String,pub parent:crate::Report,
    pub request:Request,pub steps:Vec<ActionEvidence>,pub product:ProductGraph,pub source_bases:usize,
    pub all_requested_actions_applied:bool,pub unsupported_or_unknown:bool,pub assumptions:Vec<String>,
}}
fn bad(s: &str) -> Error {
    crate::bad(&format!("processing: {s}"))
}
fn hex(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
fn call_shape(c: Call, l: &Limits) -> bool {
    c.strand < l.max_strands && c.index < l.max_source_bases
}
fn actual(s: &Substrates, c: Call) -> bool {
    s.strands
        .get(c.strand)
        .is_some_and(|st| c.index < st.sequence.len())
}
fn assumption(s: &str) -> bool {
    crate::text(s, 4096) && !s.chars().any(char::is_control)
}
pub fn request_preflight(r: &Request, l: &Limits) -> Result<(), Error> {
    crate::limits(l)?;
    let b = &r.limits;
    if b.max_actions == 0
        || b.max_actions > 64
        || b.max_product_chains == 0
        || b.max_product_chains > 128
        || b.max_mapping_segments == 0
        || b.max_mapping_segments > 16384
        || b.max_report_bytes == 0
        || b.max_report_bytes > 32 * 1024 * 1024
        || r.actions.len() > b.max_actions
        || r.actions.len() > l.max_steps
    {
        return Err(bad("processing limits/actions exceed native bounds"));
    }
    for a in &r.actions {
        match a {
            Action::CleaveThreePrimeFlap {
                bond,
                displacement_step,
                assumption: text,
                ..
            } => {
                if !call_shape(bond.upstream, l)
                    || !call_shape(bond.downstream, l)
                    || *displacement_step >= l.max_steps
                    || !assumption(text)
                {
                    return Err(bad(
                        "bounded explicit flap bond/profile/assumption required",
                    ));
                }
            }
            Action::SealNick {
                selection,
                ligase_assumption,
            } => {
                if !call_shape(selection.three, l)
                    || !call_shape(selection.five, l)
                    || !assumption(ligase_assumption)
                {
                    return Err(bad("bounded nick calls and ligase assumption required"));
                }
            }
        }
    }
    Ok(())
}
pub fn preflight(
    s: &Substrates,
    original: &[Transfer],
    l: &Limits,
    r: &Request,
) -> Result<usize, Error> {
    request_preflight(r, l)?;
    let n = crate::preflight(s, original, l)?;
    let steps = crate::add(original.len(), r.actions.len())?;
    crate::bounded(steps, l.max_steps, Budget::Steps)?;
    let cuts = r
        .actions
        .iter()
        .filter(|a| matches!(a, Action::CleaveThreePrimeFlap { .. }))
        .count();
    crate::bounded(
        crate::add(s.strands.len(), cuts)?,
        r.limits.max_product_chains,
        Budget::Strands,
    )?;
    for a in &r.actions {
        let calls = match a {
            Action::CleaveThreePrimeFlap { bond, .. } => [bond.upstream, bond.downstream],
            Action::SealNick { selection, .. } => [selection.three, selection.five],
        };
        if calls.into_iter().any(|c| !actual(s, c)) {
            return Err(bad("action call is outside an original immutable source"));
        }
    }
    let cells = crate::add(crate::mul(n, 256)?, crate::mul(steps, 256)?)?;
    crate::bounded(cells, l.max_work_cells, Budget::WorkCells)?;
    let bytes = crate::add(
        crate::mul(n, 512)?,
        crate::add(
            crate::mul(s.pairing.len(), 256)?,
            crate::add(crate::mul(steps, 16384)?, 65536)?,
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
/// Advisory extraction from a bounded report only. It grants no freshness or
/// chemical authority; simulate independently checks every actual source bond.
pub fn flap_proposals(parent: &crate::Report) -> Result<Vec<FlapProposal>, Error> {
    if parent.steps.len() > 64 {
        return Err(bad("parent transfer evidence bound"));
    }
    let mut proposals = Vec::new();
    for e in &parent.steps {
        if e.outcome != Outcome::Applied {
            continue;
        }
        if let Some(d) = &e.displacement {
            if d.kind == DisplacementKind::DuplexAnchoredThreePrimeFlap
                && d.calls_five_to_three.len() == 4
            {
                if let Some(anchor) = d.five_prime_continuation {
                    if !call_shape(anchor, &parent.limits)
                        || d.calls_five_to_three
                            .iter()
                            .any(|c| !call_shape(*c, &parent.limits))
                    {
                        return Err(bad("flap advisory call bound"));
                    }
                    proposals.push(FlapProposal {
                        displacement_step: e.step,
                        bond: Bond {
                            upstream: anchor,
                            downstream: d.calls_five_to_three[0],
                        },
                        flap_calls: d.calls_five_to_three.clone(),
                    });
                }
            }
        }
    }
    Ok(proposals)
}
fn rebuild(s: &Substrates, parent: &crate::Report, n: usize) -> Result<Graph, Error> {
    let mut g = Graph::new(s, n)?;
    for step in parent
        .steps
        .iter()
        .filter(|e| e.outcome == Outcome::Applied)
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
        let a = g.flat(step.donor_call);
        let b = g.flat(step.acceptor_call);
        g.next[a] = Some(b);
        g.prev[b] = Some(a);
    }
    Ok(g)
}
fn cut_eligibility(
    s: &Substrates,
    g: &Graph,
    parent: &crate::Report,
    bond: &Bond,
    step: usize,
) -> Result<Vec<Call>, (Outcome, &'static str)> {
    let Some(e) = parent.steps.get(step) else {
        return Err((
            Outcome::Unsupported,
            "the selected original displacement step is absent",
        ));
    };
    let Some(d) = &e.displacement else {
        return Err((
            Outcome::Unsupported,
            "selected step did not record a flap displacement",
        ));
    };
    if e.outcome != Outcome::Applied
        || d.kind != DisplacementKind::DuplexAnchoredThreePrimeFlap
        || d.calls_five_to_three.len() != 4
        || d.five_prime_continuation != Some(bond.upstream)
        || d.calls_five_to_three.first() != Some(&bond.downstream)
    {
        return Err((
            Outcome::Unsupported,
            "cut must identify the exact recorded duplex-anchored terminal four-call flap junction",
        ));
    }
    let anchor = g.flat(bond.upstream);
    let first = g.flat(bond.downstream);
    if g.next[anchor] != Some(first) || g.prev[first] != Some(anchor) {
        return Err((
            Outcome::BlockedUnderDeclaredModel,
            "the selected flap backbone bond is not currently intact or was already cut",
        ));
    }
    if g.pairs[anchor].is_none() || g.pairs[anchor].map(|x| g.call(x)) != d.continuation_pair {
        return Err((
            Outcome::Unsupported,
            "recorded flap anchor is not currently attached to its exact duplex template",
        ));
    }
    if !g.fully_concrete(s) {
        return Err((
            Outcome::Unknown,
            "ambiguous source calls do not establish this exact flap profile",
        ));
    }
    for (i, c) in d.calls_five_to_three.iter().enumerate() {
        let node = g.flat(*c);
        if g.pairs[node].is_some()
            || g.displaced[node] != Some(step)
            || g.prev[node]
                != Some(if i == 0 {
                    anchor
                } else {
                    g.flat(d.calls_five_to_three[i - 1])
                })
            || g.next[node] != (d.calls_five_to_three.get(i + 1).map(|c| g.flat(*c)))
        {
            return Err((
                Outcome::Unsupported,
                "recorded four-call flap is not an intact unpaired terminal continuation in the current state",
            ));
        }
    }
    Ok(d.calls_five_to_three.clone())
}
/// Original transfers run once; ordered actions mutate only a private reconstructed
/// state. No external graph or report can supply incoming chemical authority.
pub fn simulate(
    s: &Substrates,
    original: &[Transfer],
    l: &Limits,
    r: &Request,
) -> Result<Report, Error> {
    let fresh = crate::simulate_with_checkpoints(s, original, l, &[])?;
    process_fresh(s, original, l, r, fresh)
}
/// Only a private-constructed fresh native execution can enter this phase.
pub fn process_fresh(
    s: &Substrates,
    original: &[Transfer],
    l: &Limits,
    r: &Request,
    mut fresh: crate::Executed,
) -> Result<Report, Error> {
    let n = preflight(s, original, l, r)?;
    if fresh.report().batch_revision != crate::batch_revision(s, original, l)? {
        return Err(bad("fresh parent source/transcript differs"));
    }
    fresh.reserve_work(
        crate::add(
            crate::mul(n, 256)?,
            crate::mul(original.len() + r.actions.len(), 256)?,
        )?,
        crate::mul(n, 512)?,
    )?;
    let parent = fresh.into_report();
    let context_revision = context(&parent.batch_revision, l, r)?;
    let mut g = rebuild(s, &parent, n)?;
    let mut ends = TerminalOverrides::default();
    let mut steps = Vec::with_capacity(r.actions.len());
    for (step, a) in r.actions.iter().enumerate() {
        match a {
            Action::CleaveThreePrimeFlap {
                bond,
                displacement_step,
                profile,
                ..
            } => {
                let anchor_template = g.pairs[g.flat(bond.upstream)].map(|x| g.call(x));
                let result = cut_eligibility(s, &g, &parent, bond, *displacement_step);
                let (outcome, reason, flap_calls) = match result {
                    Ok(calls) => {
                        let u = g.flat(bond.upstream);
                        let d = g.flat(bond.downstream);
                        g.next[u] = None;
                        g.prev[d] = None;
                        ends.three.insert(u, ThreePrime::Hydroxyl);
                        ends.five.insert(d, FivePrime::Phosphate);
                        (
                            Outcome::Applied,
                            "one explicitly declared ideal hydrolytic flap-junction cut; excised source calls and original distal chemistry retained",
                            calls,
                        )
                    }
                    Err((outcome, reason)) => (outcome, reason, Vec::new()),
                };
                steps.push(ActionEvidence::FlapCleavage(FlapCleavageEvidence {
                    step,
                    bond: bond.clone(),
                    displacement_step: *displacement_step,
                    profile: *profile,
                    outcome,
                    reason: reason.into(),
                    flap_calls,
                    anchor_template,
                    created_three_prime: (outcome == Outcome::Applied)
                        .then_some(ThreePrime::Hydroxyl),
                    created_five_prime: (outcome == Outcome::Applied)
                        .then_some(FivePrime::Phosphate),
                }));
            }
            Action::SealNick { selection, .. } => {
                let u = g.flat(selection.three);
                let d = g.flat(selection.five);
                let three_prime = crate::effective_three(s, &g, u, &ends);
                let five_prime = crate::effective_five(s, &g, d, &ends);
                let template_three = g.pairs[u].map(|x| g.call(x));
                let template_five = g.pairs[d].map(|x| g.call(x));
                let (outcome, reason) = match nick::eligible_with_chemistry(
                    s,
                    &g,
                    selection,
                    three_prime.clone(),
                    five_prime,
                ) {
                    Ok(()) => {
                        g.next[u] = Some(d);
                        g.prev[d] = Some(u);
                        (
                            Outcome::Applied,
                            "one declared templated3′OH/5′P nick sealed; no pairing change, phosphorylation or implicit repair",
                        )
                    }
                    Err(e) => e,
                };
                steps.push(ActionEvidence::NickSealing(nick::SealEvidence {
                    step,
                    selection: selection.clone(),
                    outcome,
                    reason: reason.into(),
                    three_prime,
                    five_prime,
                    template_three,
                    template_five,
                }));
            }
        }
    }
    let product =
        crate::output_with_terminals(s, &g, parent.product.released_activations.clone(), &ends)?;
    let all_requested_actions_applied =
        !steps.is_empty() && steps.iter().all(|e| e.outcome() == Outcome::Applied);
    let unsupported_or_unknown = parent.unsupported_or_unknown
        || steps
            .iter()
            .any(|e| matches!(e.outcome(), Outcome::Unknown | Outcome::Unsupported));
    let report=Report{engine:ENGINE.into(),model:MODEL.into(),context_revision,parent,request:r.clone(),steps,product,source_bases:n,all_requested_actions_applied,unsupported_or_unknown,
        assumptions:vec!["The ideal hydrolytic flap profile and ligase availability/conditions are explicit operator hypotheses, not named-nuclease, commercial-kit, host-repair, kinetics or experimental-success guarantees.".into(),"Only the exact recorded duplex-anchored terminal four-call displacement is eligible for the declared cut. It creates3′OH/5′P at that bond; no source nucleotide, distal terminal chemistry, existing activation, base pair or excised material is deleted.".into(),"Only exposed3′OH/5′P calls on a continuous exact antiparallel template may be sealed. No phosphorylation, polymerase extension, dephosphorylation or implicit flap processing occurs.".into(),"All original calls partition exactly once across the complete final graph, including free excised components. Selecting one paired component is not isolation or a mixture winner. Unsupported/unknown/blocked actions retain their incoming state; subsequent actions describe that explicit branch.".into()]};
    validate_report(s, original, l, r, &report)?;
    Ok(report)
}
fn three_shape(t: &ThreePrime) -> Result<(), Error> {
    if let ThreePrime::Phosphotyrosyl { activation_id, .. } = t {
        if !crate::text(activation_id, 256) {
            return Err(bad("activation text bound"));
        }
    }
    Ok(())
}
fn product_shape(
    p: &ProductGraph,
    l: &Limits,
    b: &ProcessingLimits,
    n: usize,
) -> Result<(), Error> {
    crate::bounded(p.chains.len(), b.max_product_chains, Budget::Strands)?;
    crate::bounded(p.source_pairing.len(), l.max_pair_spans, Budget::PairSpans)?;
    crate::bounded(p.nicks.len(), b.max_mapping_segments, Budget::Strands)?;
    crate::bounded(p.unpaired.len(), n, Budget::SourceBases)?;
    crate::bounded(p.components.len(), b.max_product_chains, Budget::Strands)?;
    let mut total = 0usize;
    let mut maps = 0usize;
    for c in &p.chains {
        total = crate::add(total, c.sequence.len())?;
        crate::bounded(total, n, Budget::SourceBases)?;
        if c.sequence.alphabet() != collagene_model::Alphabet::Dna {
            return Err(bad("child DNA alphabet changed"));
        }
        maps = crate::add(maps, c.sources.len())?;
        crate::bounded(maps, b.max_mapping_segments, Budget::Strands)?;
        for m in &c.sources {
            if m.source_strand >= l.max_strands
                || m.start >= m.end
                || m.end > n
                || m.output_start
                    .checked_add(m.end - m.start)
                    .is_none_or(|e| e > c.sequence.len())
            {
                return Err(bad("source interval shape changed"));
            }
        }
        if let ChainEnds::Linear { three_prime, .. } = &c.ends {
            three_shape(three_prime)?;
        }
    }
    if total != n {
        return Err(bad("child source call count changed"));
    }
    for span in &p.source_pairing {
        if span.strand_a >= l.max_strands
            || span.strand_b >= l.max_strands
            || span.length == 0
            || span.start_a.checked_add(span.length).is_none_or(|x| x > n)
            || span.start_b.checked_add(span.length).is_none_or(|x| x > n)
        {
            return Err(bad("pair span shape changed"));
        }
    }
    for k in &p.nicks {
        if !call_shape(k.three, l) || !call_shape(k.five, l) {
            return Err(bad("nick call bounds"));
        }
        three_shape(&k.three_prime)?;
    }
    for u in &p.unpaired {
        if u.chain >= p.chains.len()
            || u.start >= u.end
            || u.end > p.chains[u.chain].sequence.len()
            || u.displaced_by.len() > l.max_steps
            || u.displaced_by.iter().any(|i| *i >= l.max_steps)
        {
            return Err(bad("retained unpaired run bounds"));
        }
    }
    for c in &p.components {
        if c.chains.len() > b.max_product_chains
            || c.chains.iter().any(|i| *i >= p.chains.len())
            || c.paired_nucleotides
                .checked_add(c.unpaired_nucleotides)
                .is_none_or(|v| v > n)
        {
            return Err(bad("component shape changed"));
        }
    }
    if p.released_activations.len() > l.max_steps
        || p.released_activations.iter().any(|s| !crate::text(s, 256))
    {
        return Err(bad("released activation bounds"));
    }
    Ok(())
}
impl Report {
    /// Readonly documentary binding/shape, never fresh chemistry authority.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        let l = &self.parent.limits;
        request_preflight(&self.request, l)?;
        crate::bounded(
            self.parent.steps.len() + self.steps.len(),
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
            || self.source_bases == 0
            || self.parent.source_revisions.is_empty()
            || self.source_bases != self.parent.source_bases
            || self.source_bases > l.max_source_bases
            || self.parent.source_revisions.len() > l.max_strands
            || self.parent.source_revisions.iter().any(|s| !hex(s))
            || self.steps.len() != self.request.actions.len()
            || self.assumptions.len() > 16
            || self.parent.assumptions.len() > 16
            || self
                .assumptions
                .iter()
                .chain(&self.parent.assumptions)
                .any(|s| s.len() > 4096)
        {
            return Err(bad("processing source/model/request binding changed"));
        }
        for (i, e) in self.parent.steps.iter().enumerate() {
            if e.step != i
                || e.reason.len() > 4096
                || e.removed_pairs.len() > 64
                || e.added_pairs.len() > 64
                || !call_shape(e.donor_call, l)
                || !call_shape(e.acceptor_call, l)
            {
                return Err(bad("parent transfer shape changed"));
            }
            for p in e.removed_pairs.iter().chain(&e.added_pairs) {
                if !call_shape(p.a, l) || !call_shape(p.b, l) {
                    return Err(bad("parent pair change shape"));
                }
            }
            if e.released_activation
                .as_ref()
                .is_some_and(|s| !crate::text(s, 256))
            {
                return Err(bad("parent activation shape"));
            }
            if let Some(d) = &e.displacement {
                if d.calls_five_to_three.len() > 64
                    || d.calls_five_to_three.iter().any(|c| !call_shape(*c, l))
                    || d.five_prime_continuation.is_some_and(|c| !call_shape(c, l))
                    || d.continuation_pair.is_some_and(|c| !call_shape(c, l))
                {
                    return Err(bad("parent displacement shape"));
                }
            }
        }
        for (i, (a, e)) in self.request.actions.iter().zip(&self.steps).enumerate() {
            if e.step() != i {
                return Err(bad("action evidence order changed"));
            }
            match (a, e) {
                (
                    Action::CleaveThreePrimeFlap {
                        bond,
                        displacement_step,
                        profile,
                        ..
                    },
                    ActionEvidence::FlapCleavage(e),
                ) => {
                    if &e.bond != bond
                        || e.displacement_step != *displacement_step
                        || e.profile != *profile
                        || e.reason.len() > 4096
                        || e.flap_calls.len() > 4
                        || e.flap_calls.iter().any(|c| !call_shape(*c, l))
                        || e.anchor_template.is_some_and(|c| !call_shape(c, l))
                        || if e.outcome == Outcome::Applied {
                            e.flap_calls.len() != 4
                                || e.created_three_prime != Some(ThreePrime::Hydroxyl)
                                || e.created_five_prime != Some(FivePrime::Phosphate)
                        } else {
                            !e.flap_calls.is_empty()
                                || e.created_three_prime.is_some()
                                || e.created_five_prime.is_some()
                        }
                    {
                        return Err(bad("flap action identity/created termini changed"));
                    }
                }
                (Action::SealNick { selection, .. }, ActionEvidence::NickSealing(e)) => {
                    if &e.selection != selection
                        || e.reason.len() > 4096
                        || e.template_three.is_some_and(|c| !call_shape(c, l))
                        || e.template_five.is_some_and(|c| !call_shape(c, l))
                    {
                        return Err(bad("nick action identity changed"));
                    }
                    three_shape(&e.three_prime)?;
                    if e.outcome == Outcome::Applied
                        && (e.three_prime != ThreePrime::Hydroxyl
                            || e.five_prime != FivePrime::Phosphate
                            || e.template_three.is_none()
                            || e.template_five.is_none())
                    {
                        return Err(bad("applied nick has unsupported chemistry"));
                    }
                }
                _ => return Err(bad("action evidence type changed")),
            }
        }
        if self.all_requested_actions_applied
            != (!self.steps.is_empty()
                && self.steps.iter().all(|e| e.outcome() == Outcome::Applied))
            || self.unsupported_or_unknown
                != (self.parent.unsupported_or_unknown
                    || self
                        .steps
                        .iter()
                        .any(|e| matches!(e.outcome(), Outcome::Unknown | Outcome::Unsupported)))
            || self.product.source_pairing != self.parent.product.source_pairing
            || self.product.released_activations != self.parent.product.released_activations
        {
            return Err(bad("processing altered pairing/activation/outcome flags"));
        }
        product_shape(
            &self.parent.product,
            l,
            &self.request.limits,
            self.source_bases,
        )?;
        product_shape(&self.product, l, &self.request.limits, self.source_bases)?;
        crate::report_size(
            self,
            l.max_report_bytes.min(self.request.limits.max_report_bytes),
        )?;
        Ok(())
    }
}
/// Every original nucleotide occurs exactly once across all final chains. Unlike
/// legacy whole-strand audits, this accepts disjoint pieces of a split source.
pub fn audit_source_partition(
    s: &Substrates,
    product: &ProductGraph,
    l: &Limits,
    b: &ProcessingLimits,
) -> Result<(), Error> {
    let n = crate::preflight(s, &[], l)?;
    product_shape(product, l, b, n)?;
    let mut offsets = Vec::with_capacity(s.strands.len());
    let mut cursor = 0usize;
    for st in &s.strands {
        offsets.push(cursor);
        cursor = crate::add(cursor, st.sequence.len())?;
    }
    let mut seen = vec![false; n];
    let mut total = 0usize;
    for chain in &product.chains {
        let mut output = 0usize;
        for m in &chain.sources {
            let source = s
                .strands
                .get(m.source_strand)
                .ok_or_else(|| bad("unknown source interval"))?;
            if m.start >= m.end || m.end > source.sequence.len() || m.output_start != output {
                return Err(bad("split source interval or output offset invalid"));
            }
            let next = crate::add(output, m.end - m.start)?;
            if chain.sequence.as_str().as_bytes().get(output..next)
                != source.sequence.as_str().as_bytes().get(m.start..m.end)
            {
                return Err(bad("processed bases differ from immutable source calls"));
            }
            for i in m.start..m.end {
                let index = offsets[m.source_strand] + i;
                if seen[index] {
                    return Err(bad("source call duplicated across processed chains"));
                }
                seen[index] = true;
            }
            output = next;
        }
        if output != chain.sequence.len() {
            return Err(bad("child chain lacks exact source interval coverage"));
        }
        total = crate::add(total, output)?;
    }
    if total != n || seen.iter().any(|seen| !*seen) {
        return Err(bad("source calls missing from retained processed mixture"));
    }
    Ok(())
}
/// Borrowed source/edge/terminal accounting, without rerunning action eligibility.
/// Fresh simulate/replay alone establishes the declared chemical hypotheses.
pub fn validate_report(
    s: &Substrates,
    original: &[Transfer],
    l: &Limits,
    r: &Request,
    report: &Report,
) -> Result<(), Error> {
    let n = preflight(s, original, l, r)?;
    report.validate_binding(&revision(s, original, l, r)?)?;
    if &report.request != r
        || &report.parent.limits != l
        || report.parent.batch_revision != crate::batch_revision(s, original, l)?
        || report.parent.steps.len() != original.len()
    {
        return Err(bad("source request or parent options changed"));
    }
    for (e, t) in report.parent.steps.iter().zip(original) {
        if &e.transfer != t
            || e.donor_call
                != (Call {
                    strand: t.donor,
                    index: s.strands[t.donor].sequence.len() - 1,
                })
            || e.acceptor_call
                != (Call {
                    strand: t.acceptor,
                    index: 0,
                })
        {
            return Err(bad("original transfer identity changed"));
        }
        for p in e.removed_pairs.iter().chain(&e.added_pairs) {
            if !actual(s, p.a) || !actual(s, p.b) {
                return Err(bad("parent pair call absent"));
            }
        }
        if let Some(d) = &e.displacement {
            if d.calls_five_to_three.iter().any(|c| !actual(s, *c))
                || d.five_prime_continuation.is_some_and(|c| !actual(s, c))
                || d.continuation_pair.is_some_and(|c| !actual(s, c))
            {
                return Err(bad("parent displacement call absent"));
            }
        }
    }
    crate::audit_sources(s, &report.parent)?;
    let mut g = rebuild(s, &report.parent, n)?;
    let mut ends = TerminalOverrides::default();
    if crate::output(s, &g, report.parent.product.released_activations.clone())?
        != report.parent.product
    {
        return Err(bad(
            "parent graph disagrees with its declared applied transfer edges",
        ));
    }
    for e in &report.steps {
        match e {
            ActionEvidence::FlapCleavage(e) => {
                if e.outcome != Outcome::Applied {
                    continue;
                }
                let d = report
                    .parent
                    .steps
                    .get(e.displacement_step)
                    .and_then(|e| e.displacement.as_ref())
                    .ok_or_else(|| bad("applied cut has no original displacement"))?;
                if d.kind != DisplacementKind::DuplexAnchoredThreePrimeFlap
                    || d.calls_five_to_three != e.flap_calls
                    || d.five_prime_continuation != Some(e.bond.upstream)
                    || e.flap_calls.first() != Some(&e.bond.downstream)
                {
                    return Err(bad("applied cut changed original flap identity"));
                }
                let u = g.flat(e.bond.upstream);
                let down = g.flat(e.bond.downstream);
                if g.next[u] != Some(down)
                    || g.prev[down] != Some(u)
                    || e.anchor_template != g.pairs[u].map(|x| g.call(x))
                {
                    return Err(bad(
                        "applied cut removed an absent bond or changed its template",
                    ));
                }
                g.next[u] = None;
                g.prev[down] = None;
                ends.three.insert(u, ThreePrime::Hydroxyl);
                ends.five.insert(down, FivePrime::Phosphate);
            }
            ActionEvidence::NickSealing(e) => {
                let u = g.flat(e.selection.three);
                let down = g.flat(e.selection.five);
                if e.three_prime != crate::effective_three(s, &g, u, &ends)
                    || e.five_prime != crate::effective_five(s, &g, down, &ends)
                    || e.template_three != g.pairs[u].map(|x| g.call(x))
                    || e.template_five != g.pairs[down].map(|x| g.call(x))
                {
                    return Err(bad(
                        "nick evidence changed actual endpoint/template chemistry",
                    ));
                }
                if e.outcome == Outcome::Applied {
                    if g.next[u].is_some() || g.prev[down].is_some() {
                        return Err(bad("applied seal reused a consumed terminus"));
                    }
                    g.next[u] = Some(down);
                    g.prev[down] = Some(u);
                }
            }
        }
    }
    if crate::output_with_terminals(
        s,
        &g,
        report.parent.product.released_activations.clone(),
        &ends,
    )? != report.product
    {
        return Err(bad(
            "processed graph differs from recorded bond and terminal changes",
        ));
    }
    audit_source_partition(s, &report.product, l, &r.limits)?;
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
        return Err(bad("report differs from fresh ordered chemical processing"));
    }
    Ok(())
}
