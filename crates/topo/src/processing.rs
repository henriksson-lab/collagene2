//! Explicit ideal flap cleavage and nick sealing inside ONE fresh private graph cursor.
//! Cleavage is an operator-declared hydrolytic hypothesis, not a nuclease/kit model.
use crate::{
    Call, DisplacementKind, Error, FivePrime, Graph, Limits, Outcome, Substrates,
    TerminalOverrides, ThreePrime, nick_sealing as nick,
};
use serde::{Deserialize, Serialize};
pub const ENGINE: &str = "topo-source-processing-2";
pub const MODEL: &str = "ideal-anchored-three-flap-hydrolysis-and-templated-nick-sealing-1";
macro_rules! dto {($($t:tt)*)=>{#[derive(Debug,Clone,PartialEq,Eq,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct Bond{pub upstream:Call,pub downstream:Call}}
dto! {#[derive(Copy)] pub enum FlapCleavageProfile{IdealDuplexAnchoredThreePrimeHydrolysisV1}}
dto! {pub enum Action{
    CleaveThreePrimeFlap{bond:Bond,displacement_step:usize,profile:FlapCleavageProfile,assumption:String},
    SealNick{selection:nick::NickSelection,ligase_assumption:String},
}}
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

impl Action {
    pub fn resource_preflight(&self) -> Result<(), Error> {
        let hypothesis = match self {
            Self::CleaveThreePrimeFlap { assumption, .. } => assumption,
            Self::SealNick {
                ligase_assumption, ..
            } => ligase_assumption,
        };
        if hypothesis.len() > 4096 {
            return Err(crate::bad("action hypothesis exceeds resource bound"));
        }
        Ok(())
    }
}
pub fn action_preflight(action: &Action, l: &Limits) -> Result<(), Error> {
    action.resource_preflight()?;
    let call = |c: Call| c.strand < l.max_strands && c.index < l.max_source_bases;
    let hypothesis = |text: &str| !text.trim().is_empty() && !text.chars().any(char::is_control);
    let valid = match action {
        Action::CleaveThreePrimeFlap {
            bond,
            displacement_step,
            assumption,
            ..
        } => {
            call(bond.upstream)
                && call(bond.downstream)
                && *displacement_step < l.max_steps
                && hypothesis(assumption)
        }
        Action::SealNick {
            selection,
            ligase_assumption,
        } => call(selection.three) && call(selection.five) && hypothesis(ligase_assumption),
    };
    if !valid {
        return Err(crate::bad(
            "bounded explicit action/call/conditions required",
        ));
    }
    Ok(())
}
pub(crate) fn calls(action: &Action) -> [Call; 2] {
    match action {
        Action::CleaveThreePrimeFlap { bond, .. } => [bond.upstream, bond.downstream],
        Action::SealNick { selection, .. } => [selection.three, selection.five],
    }
}
/// Advisory evidence only. Native execution independently checks the exact bond.
pub fn flap_proposals(parent: &crate::Report) -> Result<Vec<FlapProposal>, Error> {
    if parent.steps.len() > 64 {
        return Err(crate::bad("bounded reaction evidence required"));
    }
    let mut proposals = Vec::new();
    for step in &parent.steps {
        if let crate::StepEvidence::Transfer(e) = step {
            if e.outcome == Outcome::Applied {
                if let Some(d) = &e.displacement {
                    if d.kind == DisplacementKind::DuplexAnchoredThreePrimeFlap
                        && d.calls_five_to_three.len() == 4
                    {
                        if let Some(anchor) = d.five_prime_continuation {
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
        }
    }
    Ok(proposals)
}
pub(crate) fn apply(
    s: &Substrates,
    g: &mut Graph,
    ends: &mut TerminalOverrides,
    prior: &[crate::StepEvidence],
    step: usize,
    action: &Action,
) -> ActionEvidence {
    match action {
        Action::CleaveThreePrimeFlap {
            bond,
            displacement_step,
            profile,
            ..
        } => {
            let anchor_template = g.pairs[g.flat(bond.upstream)].map(|x| g.call(x));
            let result = cut_eligibility(s, g, prior, bond, *displacement_step);
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
            ActionEvidence::FlapCleavage(FlapCleavageEvidence {
                step,
                bond: bond.clone(),
                displacement_step: *displacement_step,
                profile: *profile,
                outcome,
                reason: reason.into(),
                flap_calls,
                anchor_template,
                created_three_prime: (outcome == Outcome::Applied).then_some(ThreePrime::Hydroxyl),
                created_five_prime: (outcome == Outcome::Applied).then_some(FivePrime::Phosphate),
            })
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
            ActionEvidence::NickSealing(nick::SealEvidence {
                step,
                selection: selection.clone(),
                outcome,
                reason: reason.into(),
                three_prime,
                five_prime,
                template_three,
                template_five,
            })
        }
    }
}
fn cut_eligibility(
    s: &Substrates,
    g: &Graph,
    prior: &[crate::StepEvidence],
    bond: &Bond,
    step: usize,
) -> Result<Vec<Call>, (Outcome, &'static str)> {
    let Some(crate::StepEvidence::Transfer(e)) = prior.get(step) else {
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
