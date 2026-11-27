//! Exact continuous-template 3′OH/5′P nick eligibility within the canonical cursor.
//! No phosphorylation, polymerase, activation consumption or host repair is inferred.
use crate::{Call, FivePrime, Graph, Outcome, Substrates, ThreePrime};
use serde::{Deserialize, Serialize};
macro_rules! dto{($($t:tt)*)=>{#[derive(Clone,Debug,PartialEq,Eq,Serialize,Deserialize)]#[cfg_attr(feature="json-schema",derive(schemars::JsonSchema))]#[serde(deny_unknown_fields)]$($t)*};}
dto! {pub struct NickSelection{pub three:Call,pub five:Call}}
dto! {pub struct SealEvidence{pub step:usize,pub selection:NickSelection,pub outcome:Outcome,pub reason:String,
pub three_prime:ThreePrime,pub five_prime:FivePrime,pub template_three:Option<Call>,pub template_five:Option<Call>}}
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
