//! Explicit element roles and directed source geometry; no transcription model.
use super::*;

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
pub enum ElementRole {
    Promoter,
    Terminator,
    PolyadenylationSignal,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ChainElement {
    pub feature_index: usize,
    pub reference_id: Option<String>,
    pub role: ElementRole,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ElementChain {
    pub id: String,
    pub cds_feature: usize,
    pub upstream: Option<ChainElement>,
    pub downstream: Option<ChainElement>,
    /// Optional explicit user constraints, not invented biological thresholds.
    pub max_upstream_distance: Option<usize>,
    pub max_downstream_distance: Option<usize>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ChainSideEvidence {
    pub feature_index: Option<usize>,
    pub role: Option<ElementRole>,
    pub reference_id: Option<String>,
    pub verdict: Verdict,
    pub distance_bases: Option<usize>,
    pub reverse: Option<bool>,
    pub reason: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "json-schema", derive(schemars::JsonSchema))]
#[serde(deny_unknown_fields)]
pub struct ElementChainReport {
    pub id: String,
    pub cds_feature: usize,
    pub upstream: ChainSideEvidence,
    pub downstream: ChainSideEvidence,
    /// Conjunction of this declaration's structural predicates, not expression.
    pub verdict: Verdict,
    pub reason: String,
}

impl ElementChain {
    pub(super) fn validate(&self, m: &Molecule) -> Result<(), Error> {
        if !text(&self.id, 128)
            || self.cds_feature >= m.features.len()
            || self.max_upstream_distance.is_some_and(|n| n > MAX_BASES)
            || self.max_downstream_distance.is_some_and(|n| n > MAX_BASES)
        {
            return Err(bad("invalid element-chain declaration"));
        }
        for e in [&self.upstream, &self.downstream].into_iter().flatten() {
            if e.feature_index >= m.features.len()
                || e.reference_id.as_ref().is_some_and(|id| !text(id, 128))
            {
                return Err(bad("invalid declared element feature/reference"));
            }
        }
        if self
            .upstream
            .as_ref()
            .is_some_and(|e| e.role != ElementRole::Promoter)
            || self
                .downstream
                .as_ref()
                .is_some_and(|e| e.role == ElementRole::Promoter)
        {
            return Err(bad(
                "upstream role must be promoter; downstream must be terminator/polyadenylation signal",
            ));
        }
        Ok(())
    }
}

fn side(
    checks: &Checks<'_>,
    chain: &ElementChain,
    element: Option<&ChainElement>,
    upstream: bool,
) -> Result<ChainSideEvidence, Error> {
    let mut result = ChainSideEvidence {
        feature_index: element.map(|e| e.feature_index),
        role: element.map(|e| e.role.clone()),
        reference_id: element.and_then(|e| e.reference_id.clone()),
        verdict: Verdict::Indeterminate,
        distance_bases: None,
        reverse: None,
        reason: String::new(),
    };
    let Some(element) = element else {
        result.reason = "No element is declared on this side; missing annotation does not establish biological absence.".into();
        return Ok(result);
    };
    if !checks.assumed(chain.cds_feature) || !checks.assumed(element.feature_index) {
        result.reason = "Both CDS and element roles need current source/annotation-pinned assumptions; unreviewed or undeclared annotation cannot establish the chain.".into();
        return Ok(result);
    }
    if !checks.m.features[chain.cds_feature]
        .kind
        .eq_ignore_ascii_case("CDS")
    {
        result.reason = "The declared coding feature is not a CDS annotation.".into();
        return Ok(result);
    }
    let Some(reference_id) = &element.reference_id else {
        result.reason = "The element has no pinned literal reference; its label is not sequence identity evidence.".into();
        return Ok(result);
    };
    let identity = geometry::reference_identity(checks, element.feature_index, reference_id)?;
    if identity.verdict != Verdict::Pass {
        // A mismatching reference invalidates the role premise, rather than
        // proving that this construct has no promoter or termination element.
        result.reason = format!(
            "Declared element identity is unresolved: {}",
            identity.reason
        );
        return Ok(result);
    }
    let Some(((cs, ce, cr), (es, ee, er))) =
        geometry::single(&checks.m.features[chain.cds_feature].location, false).zip(
            geometry::single(&checks.m.features[element.feature_index].location, false),
        )
    else {
        result.reason = "Compound, fuzzy or unresolved coordinate paths do not define the simple directed interval required by this chain predicate.".into();
        return Ok(result);
    };
    result.reverse = Some(er);
    if cr != er {
        result.verdict = Verdict::Fail;
        result.reason = "The exact declared element and CDS are on opposite strands; this is an orientation fact, not a functional verdict.".into();
        return Ok(result);
    }
    if es < ce && cs < ee {
        result.verdict = Verdict::Fail;
        result.reason = "Declared element and CDS intervals overlap, violating this explicit flanking-chain requirement.".into();
        return Ok(result);
    }
    let from_left = upstream != cr;
    let distance = if checks.m.topology == Topology::Circular {
        let n = checks.m.sequence.len();
        Some(if from_left {
            (cs + n - ee) % n
        } else {
            (es + n - ce) % n
        })
    } else if from_left {
        cs.checked_sub(ee)
    } else {
        es.checked_sub(ce)
    };
    let Some(distance) = distance else {
        result.verdict = Verdict::Fail;
        result.reason =
            "The element lies on the wrong side of the declared CDS on this linear source.".into();
        return Ok(result);
    };
    result.distance_bases = Some(distance);
    let maximum = if upstream {
        chain.max_upstream_distance
    } else {
        chain.max_downstream_distance
    };
    result.verdict = if maximum.is_some_and(|n| distance > n) {
        Verdict::Fail
    } else {
        Verdict::Pass
    };
    result.reason = format!(
        "Same-strand directed interval distance is {distance} bases{}. This proves declared order/orientation/spacing only; the element role and biological activity are not verified.",
        maximum
            .map(|n| format!("; user-declared maximum is {n}"))
            .unwrap_or_default(),
    );
    Ok(result)
}

pub(super) fn inspect(checks: &mut Checks<'_>) -> Result<Vec<ElementChainReport>, Error> {
    let mut reports = Vec::new();
    let m = checks.m;
    let o = checks.o;
    for (i, f) in m.features.iter().enumerate() {
        if f.kind.eq_ignore_ascii_case("CDS")
            && !o.element_chains.iter().any(|c| c.cds_feature == i)
        {
            checks.push(Rule::ElementChain, Verdict::Indeterminate, &[i],
                "No explicit flanking-element chain is declared for this CDS. Missing declarations/annotations do not establish promoter or termination-element absence.")?;
        }
    }
    for chain in &o.element_chains {
        let upstream = side(checks, chain, chain.upstream.as_ref(), true)?;
        let downstream = side(checks, chain, chain.downstream.as_ref(), false)?;
        for (rule, evidence) in [
            (Rule::UpstreamElement, &upstream),
            (Rule::DownstreamElement, &downstream),
        ] {
            let mut indices = vec![chain.cds_feature];
            if let Some(i) = evidence.feature_index {
                indices.push(i);
            }
            checks.push(rule, evidence.verdict.clone(), &indices, &evidence.reason)?;
        }
        let mut verdict =
            if upstream.verdict == Verdict::Fail || downstream.verdict == Verdict::Fail {
                Verdict::Fail
            } else if upstream.verdict == Verdict::Indeterminate
                || downstream.verdict == Verdict::Indeterminate
            {
                Verdict::Indeterminate
            } else {
                Verdict::Pass
            };
        let mut reason = "This verdict concerns this declared flanking chain only; no promoter, processing, termination or expression activity is inferred.".to_string();
        if verdict == Verdict::Pass && checks.m.topology == Topology::Circular {
            let c = geometry::single(&checks.m.features[chain.cds_feature].location, false)
                .ok_or_else(|| bad("missing resolved CDS geometry"))?;
            let lengths = [&chain.upstream, &chain.downstream]
                .into_iter()
                .flatten()
                .try_fold(0usize, |n, e| {
                    let (s, e, _) =
                        geometry::single(&checks.m.features[e.feature_index].location, false)
                            .ok_or_else(|| bad("missing resolved flanking-element geometry"))?;
                    n.checked_add(e - s).ok_or(Error::TooLarge)
                })?;
            let span = lengths + c.1 - c.0
                + upstream.distance_bases.unwrap_or(0)
                + downstream.distance_bases.unwrap_or(0);
            if span > checks.m.sequence.len() {
                verdict = Verdict::Fail;
                reason = "Directed upstream→CDS→downstream chain exceeds one source turn or reuses overlapping flanks; the declared circular order is inconsistent.".into();
            }
        }
        let mut indices = vec![chain.cds_feature];
        indices.extend(
            [upstream.feature_index, downstream.feature_index]
                .into_iter()
                .flatten(),
        );
        checks.push(Rule::ElementChain, verdict.clone(), &indices, &reason)?;
        reports.push(ElementChainReport {
            id: chain.id.clone(),
            cds_feature: chain.cds_feature,
            upstream,
            downstream,
            verdict,
            reason,
        });
    }
    Ok(reports)
}

pub(super) fn validate_reports(report: &Report) -> Result<(), Error> {
    let mut ids = BTreeSet::new();
    for chain in &report.element_chains {
        if !text(&chain.id, 128)
            || !ids.insert(&chain.id)
            || chain.cds_feature >= report.source_feature_count
            || !text(&chain.reason, 4096)
        {
            return Err(bad("invalid retained element-chain binding"));
        }
        for side in [&chain.upstream, &chain.downstream] {
            if !text(&side.reason, 4096)
                || side
                    .feature_index
                    .is_some_and(|i| i >= report.source_feature_count)
                || side.reference_id.as_ref().is_some_and(|s| !text(s, 128))
                || side
                    .distance_bases
                    .is_some_and(|n| n > report.source_length)
                || (side.feature_index.is_none()
                    && (side.role.is_some()
                        || side.reference_id.is_some()
                        || side.distance_bases.is_some()
                        || side.reverse.is_some()
                        || side.verdict != Verdict::Indeterminate))
                || (side.verdict == Verdict::Pass
                    && (side.reference_id.is_none()
                        || side.distance_bases.is_none()
                        || side.reverse.is_none()
                        || side.role.is_none()
                        || !report
                            .references
                            .iter()
                            .any(|r| Some(&r.id) == side.reference_id.as_ref())))
            {
                return Err(bad("invalid retained directed-element geometry"));
            }
        }
        if chain.verdict == Verdict::Pass
            && (chain.upstream.verdict != Verdict::Pass
                || chain.downstream.verdict != Verdict::Pass)
        {
            return Err(bad("passing declared chain lacks passing flanks"));
        }
    }
    Ok(())
}
