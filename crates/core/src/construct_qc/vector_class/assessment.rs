//! Data-driven checks of independently bound role premises and directed geometry.
use super::*;
use profiles::Orientation;

data! {
    #[derive(Copy)]
    pub enum Rule {
        ReferenceIntegrity,
        RolePremise,
        DeclaredMultiplicity,
        BoundarySelection,
        PayloadMeasurement,
        Containment,
        Orientation,
        Order,
        DeclaredSpacing,
        EmpiricalContext,
    }
}
data! {
    pub struct Span {
        pub start: usize,
        pub end: usize,
        pub reverse: bool,
    }
}
data! {
    pub struct BindingEvidence {
        pub binding: RoleBinding,
        /// The actual canonical literal-reference claim result is retained even
        /// if the role premise and its dependent geometry remain indeterminate.
        pub reference_verdict: Verdict,
        pub premise: Verdict,
        pub reason: String,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub interval: Option<Span>,
    }
}
data! {
    pub struct BoundaryMeasurement {
        pub selection: BoundarySelection,
        pub start_bond: usize,
        pub end_bond: usize,
        /// Includes both explicitly selected boundary elements.
        pub outer_bases: usize,
        /// Excludes both boundary elements; distinct from total plasmid length.
        pub interior_bases: usize,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Span, 2>")]
        pub source_spans: Vec<Span>,
    }
}
data! {
    pub struct Finding {
        pub rule: Rule,
        pub verdict: Verdict,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, String, 4>")]
        pub binding_ids: Vec<String>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, usize, 4>")]
        pub feature_indices: Vec<usize>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Span, 4>")]
        pub coordinates: Vec<Span>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, String, 8>")]
        pub source_ids: Vec<String>,
        pub message: String,
    }
}
data! {
    pub struct Report {
        pub engine_version: String,
        pub context_revision: String,
        pub source_revision: String,
        pub source_length: usize,
        pub source_feature_count: usize,
        pub source_topology: Topology,
        pub annotations_revision: String,
        pub assumptions_revision: String,
        pub references_revision: String,
        pub request: Request,
        pub profile: Profile,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, BindingEvidence, 64>")]
        pub bindings: Vec<BindingEvidence>,
        #[serde(deserialize_with = "crate::encoding::required_option")]
        #[cfg_attr(feature = "json-schema", schemars(required))]
        pub boundaries: Option<BoundaryMeasurement>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, Finding, 256>")]
        pub findings: Vec<Finding>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, super::super::required_parts::RequiredPartReport, 64>")]
        pub literal_parts: Vec<super::super::required_parts::RequiredPartReport>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, super::super::required_parts::ReferenceSearch, 64>")]
        pub literal_searches: Vec<super::super::required_parts::ReferenceSearch>,
        #[serde(deserialize_with = "collagene_model::decode::bounded_vec::<_, super::super::Finding, 64>")]
        pub literal_findings: Vec<super::super::Finding>,
        pub literal_search_work: usize,
        pub statement: String,
        pub report_revision: String,
    }
}
struct Findings(Vec<Finding>);
impl Findings {
    fn push(
        &mut self,
        rule: Rule,
        verdict: Verdict,
        nodes: &[&BindingEvidence],
        sources: &[String],
        message: impl Into<String>,
    ) -> Result<(), Error> {
        if self.0.len() >= 256 || nodes.len() > 4 || sources.len() > 8 {
            return Err(Error::TooLarge);
        }
        let message = message.into();
        if message.len() > 2048 {
            return Err(Error::TooLarge);
        }
        self.0.push(Finding {
            rule,
            verdict,
            message,
            binding_ids: nodes.iter().map(|node| node.binding.id.clone()).collect(),
            feature_indices: nodes
                .iter()
                .map(|node| node.binding.feature_index)
                .collect(),
            coordinates: nodes
                .iter()
                .filter_map(|node| node.interval.clone())
                .collect(),
            source_ids: sources.to_vec(),
        });
        Ok(())
    }
}
fn identity(
    source: &Molecule,
    request: &Request,
    checks: &super::super::Checks<'_>,
    findings: &mut Findings,
) -> Result<Vec<BindingEvidence>, Error> {
    request.bindings.iter().map(|binding| {
        let assumed = checks.assumed(binding.feature_index);
        let (reference_verdict, reason) = if assumed {
            let evidence = super::super::geometry::reference_identity(checks, binding.feature_index, &binding.reference_id)?;
            (evidence.verdict, evidence.reason)
        } else {
            (Verdict::Indeterminate, "Role requires a current source/annotation-pinned assumption; missing, stale or review-required annotation does not establish this premise.".into())
        };
        let qualified = assumed && reference_verdict == Verdict::Pass;
        let interval = if qualified {
            super::super::geometry::single(&source.features[binding.feature_index].location, false)
                .map(|(start, end, reverse)| Span { start, end, reverse })
        } else { None };
        let result = BindingEvidence {
            binding: binding.clone(), reference_verdict: reference_verdict.clone(),
            premise: if qualified { Verdict::Pass } else { Verdict::Indeterminate },
            reason: reason.clone(), interval,
        };
        findings.push(Rule::ReferenceIntegrity, reference_verdict, &[&result], &[], reason)?;
        if !qualified || result.interval.is_none() {
            findings.push(Rule::RolePremise, Verdict::Indeterminate, &[&result], &[],
                "Dependent role geometry is unresolved. A mismatching literal claim is retained as an identity failure, not interpreted as biological absence. Compound/fuzzy paths have no guessed simple interval.")?;
        }
        Ok(result)
    }).collect()
}
fn boundary(
    source: &Molecule,
    request: &Request,
    profile: &Profile,
    nodes: &[BindingEvidence],
    findings: &mut Findings,
) -> Result<Option<BoundaryMeasurement>, Error> {
    let Some(selection) = &request.boundaries else {
        findings.push(Rule::BoundarySelection, Verdict::Indeterminate, &[], &[], "No exact boundary bindings/axis were declared; payload length is not replaced by total plasmid size.")?;
        return Ok(None);
    };
    let start = nodes
        .iter()
        .find(|node| node.binding.id == selection.start_binding)
        .ok_or_else(|| bad("missing selected start binding"))?;
    let end = nodes
        .iter()
        .find(|node| node.binding.id == selection.end_binding)
        .ok_or_else(|| bad("missing selected end binding"))?;
    let Some((first, last)) = start.interval.as_ref().zip(end.interval.as_ref()) else {
        findings.push(Rule::BoundarySelection, Verdict::Indeterminate, &[start, end], &[], "Both boundary role premises and simple intervals must be established; no payload extent is inferred from labels alone.")?;
        return Ok(None);
    };
    if start.binding.role_id != profile.definition.start_role
        || end.binding.role_id != profile.definition.end_role
    {
        findings.push(
            Rule::BoundarySelection,
            Verdict::Fail,
            &[start, end],
            &[],
            "Explicit endpoint roles do not match the installed profile's declared endpoint roles.",
        )?;
        return Ok(None);
    }
    let n = source.sequence.len();
    let start_bond = if selection.reverse {
        first.end
    } else {
        first.start
    };
    let end_bond = if selection.reverse {
        last.start
    } else {
        last.end
    };
    let (start_bond, end_bond) = if source.topology == Topology::Circular {
        (start_bond % n, end_bond % n)
    } else {
        (start_bond, end_bond)
    };
    let distance = if source.topology == Topology::Circular {
        let distance = if selection.reverse {
            (start_bond + n - end_bond) % n
        } else {
            (end_bond + n - start_bond) % n
        };
        // Both distinct boundary elements must be included. Coincident outer
        // bonds therefore mean one full turn, never an empty payload span.
        Some(if distance == 0 { n } else { distance })
    } else if selection.reverse {
        start_bond.checked_sub(end_bond)
    } else {
        end_bond.checked_sub(start_bond)
    };
    let Some(outer_bases) = distance.filter(|length| *length > 0) else {
        findings.push(Rule::BoundarySelection, Verdict::Fail, &[start, end], &[], "The declared endpoints/axis do not specify a positive single-turn span on this source.")?;
        return Ok(None);
    };
    let Some(interior_bases) = outer_bases
        .checked_sub(first.end - first.start)
        .and_then(|length| length.checked_sub(last.end - last.start))
    else {
        findings.push(
            Rule::BoundarySelection,
            Verdict::Fail,
            &[start, end],
            &[],
            "The selected outer span cannot contain both nonoverlapping boundary elements.",
        )?;
        return Ok(None);
    };
    let source_spans = if !selection.reverse {
        if start_bond < end_bond {
            vec![Span {
                start: start_bond,
                end: end_bond,
                reverse: false,
            }]
        } else {
            [
                Span {
                    start: start_bond,
                    end: n,
                    reverse: false,
                },
                Span {
                    start: 0,
                    end: end_bond,
                    reverse: false,
                },
            ]
            .into_iter()
            .filter(|span| span.start < span.end)
            .collect()
        }
    } else if end_bond < start_bond {
        vec![Span {
            start: end_bond,
            end: start_bond,
            reverse: true,
        }]
    } else {
        [
            Span {
                start: 0,
                end: start_bond,
                reverse: true,
            },
            Span {
                start: end_bond,
                end: n,
                reverse: true,
            },
        ]
        .into_iter()
        .filter(|span| span.start < span.end)
        .collect()
    };
    findings.push(Rule::PayloadMeasurement, Verdict::Informational, &[start, end], &[], format!("Declared outer boundary span is {outer_bases} bases; interior excluding both boundary elements is {interior_bases} bases. Total source length {n} is not substituted. This is a coordinate measurement, not a packaging budget verdict."))?;
    Ok(Some(BoundaryMeasurement {
        selection: selection.clone(),
        start_bond,
        end_bond,
        outer_bases,
        interior_bases,
        source_spans,
    }))
}
fn placement(
    node: &BindingEvidence,
    boundary: &BoundaryMeasurement,
    source: &Molecule,
) -> Option<(usize, usize)> {
    let interval = node.interval.as_ref()?;
    let distance = if source.topology == Topology::Circular {
        Some(if boundary.selection.reverse {
            (boundary.start_bond + source.sequence.len() - interval.end) % source.sequence.len()
        } else {
            (interval.start + source.sequence.len() - boundary.start_bond) % source.sequence.len()
        })
    } else if boundary.selection.reverse {
        boundary.start_bond.checked_sub(interval.end)
    } else {
        interval.start.checked_sub(boundary.start_bond)
    }?;
    let end = distance.checked_add(interval.end - interval.start)?;
    (end <= boundary.outer_bases).then_some((distance, end))
}
fn gap(
    first: &BindingEvidence,
    last: &BindingEvidence,
    boundary: Option<&BoundaryMeasurement>,
    source: &Molecule,
) -> Option<usize> {
    let boundary = boundary?;
    let (_, first_end) = placement(first, boundary, source)?;
    let (last_start, _) = placement(last, boundary, source)?;
    last_start.checked_sub(first_end)
}
fn check_order(
    findings: &mut Findings,
    rule: Rule,
    nodes: &[&BindingEvidence],
    boundary: Option<&BoundaryMeasurement>,
    source: &Molecule,
    minimum: usize,
    maximum: Option<usize>,
    sources: &[String],
    rationale: &str,
) -> Result<(), Error> {
    if nodes.len() != 2
        || nodes
            .iter()
            .any(|node| node.premise != Verdict::Pass || node.interval.is_none())
        || boundary.is_none()
    {
        return findings.push(rule, Verdict::Indeterminate, &nodes[..nodes.len().min(4)], sources,
            format!("Directed spacing/order needs two uniquely selected qualified simple intervals and explicit boundaries. {rationale}"));
    }
    let measured = gap(nodes[0], nodes[1], boundary, source);
    findings.push(rule, if measured.is_some_and(|gap| gap >= minimum && maximum.is_none_or(|maximum| gap <= maximum)) { Verdict::Pass } else { Verdict::Fail }, nodes, sources,
        format!("Directed gap is {measured:?} bases; declared minimum {minimum}, maximum {maximum:?}. A missing directed gap here means overlap, wrong order or outside-boundary geometry, not absent biology. {rationale}"))
}
pub fn assess(source: &Molecule, request: &Request) -> Result<Report, Error> {
    let context_revision = revision(source, request)?;
    let profile = profiles::resolve(&request.profile)?;
    let options = request.options();
    // One canonical assumption/reference-facts owner, shared with ordinary QC.
    let mut checks = super::super::Checks::new(source, &options)?;
    let mut findings = Findings(vec![]);
    let nodes = identity(source, request, &checks, &mut findings)?;
    let boundaries = boundary(source, request, &profile, &nodes, &mut findings)?;
    for role in &profile.definition.roles {
        let selected = nodes
            .iter()
            .filter(|node| node.binding.role_id == role.id)
            .collect::<Vec<_>>();
        let qualified = selected
            .iter()
            .filter(|node| node.premise == Verdict::Pass)
            .count();
        let verdict = if selected.is_empty() || qualified != selected.len() {
            Verdict::Indeterminate
        } else if qualified < role.minimum || qualified > role.maximum {
            Verdict::Fail
        } else {
            Verdict::Pass
        };
        findings.push(Rule::DeclaredMultiplicity, verdict, &selected[..selected.len().min(4)], &role.source_ids,
            format!("Role '{}': {qualified} qualified of {} explicitly bound node(s), required {}..{}. Missing/unqualified annotations do not prove absence; these counts are not exhaustive biological copy numbers.", role.id, selected.len(), role.minimum, role.maximum))?;
        for node in selected {
            if let Some(boundary) = &boundaries {
                let positioned = placement(node, boundary, source);
                if role.within_boundaries {
                    findings.push(Rule::Containment, if node.interval.is_none() { Verdict::Indeterminate } else if positioned.is_some() { Verdict::Pass } else { Verdict::Fail }, &[node], &role.source_ids,
                        "Check this independently qualified declared element within the selected outer span; uncertain role premises remain unresolved.")?;
                }
                let orientation = node
                    .interval
                    .as_ref()
                    .map(|interval| match role.orientation {
                        Orientation::SameAsAxis => interval.reverse == boundary.selection.reverse,
                        Orientation::OppositeToAxis => {
                            interval.reverse != boundary.selection.reverse
                        }
                        Orientation::Either => true,
                    });
                findings.push(Rule::Orientation, match orientation { Some(true) => Verdict::Pass, Some(false) => Verdict::Fail, None => Verdict::Indeterminate }, &[node], &role.source_ids,
                    format!("Declared feature orientation is checked relative to the selected graph axis under {:?}; this does not authenticate element role or physical structure.", role.orientation))?;
            }
        }
    }
    for order in &profile.definition.order {
        let first = nodes
            .iter()
            .filter(|node| node.binding.role_id == order.upstream_role)
            .collect::<Vec<_>>();
        let last = nodes
            .iter()
            .filter(|node| node.binding.role_id == order.downstream_role)
            .collect::<Vec<_>>();
        let selected = if first.len() == 1 && last.len() == 1 {
            vec![first[0], last[0]]
        } else {
            vec![]
        };
        check_order(
            &mut findings,
            Rule::Order,
            &selected,
            boundaries.as_ref(),
            source,
            order.minimum_gap,
            order.maximum_gap,
            &order.source_ids,
            "Installed declared-graph rule; no transcription/activity prediction.",
        )?;
    }
    for spacing in &request.spacing {
        let first = nodes
            .iter()
            .find(|node| node.binding.id == spacing.upstream_binding)
            .ok_or_else(|| bad("missing upstream spacing binding"))?;
        let last = nodes
            .iter()
            .find(|node| node.binding.id == spacing.downstream_binding)
            .ok_or_else(|| bad("missing downstream spacing binding"))?;
        check_order(
            &mut findings,
            Rule::DeclaredSpacing,
            &[first, last],
            boundaries.as_ref(),
            source,
            spacing.minimum_gap,
            spacing.maximum_gap,
            &[],
            &spacing.rationale,
        )?;
    }
    for context in &profile.definition.numeric_context {
        findings.push(Rule::EmpiricalContext, Verdict::Informational, &[], &context.source_ids,
            format!("{}: reported {}..{} {}. {} {} This context is never applied as a pass/fail threshold to the current source.", context.metric, context.minimum, context.maximum, context.unit, context.population, context.interpretation))?;
    }
    let (literal_parts, literal_searches, literal_search_work, literal_findings) =
        if request.literal_parts.is_empty() {
            (vec![], vec![], 0, vec![])
        } else {
            let result = super::super::required_parts::inspect(&mut checks)?;
            (
                result.reports,
                result.searches,
                result.work,
                checks.findings,
            )
        };
    let mut report = Report {
        engine_version: ENGINE_VERSION.into(),
        context_revision,
        source_revision: source.revision(),
        source_length: source.sequence.len(),
        source_feature_count: source.features.len(),
        source_topology: source.topology,
        annotations_revision: hash(&source.features)?,
        assumptions_revision: hash(&request.annotations)?,
        references_revision: hash(&request.references)?,
        request: request.clone(),
        profile,
        bindings: nodes,
        boundaries,
        findings: findings.0,
        literal_parts,
        literal_searches,
        literal_findings,
        literal_search_work,
        statement: STATEMENT.into(),
        report_revision: String::new(),
    };
    bounded(&report, MAX_REPORT_BYTES)?;
    report.report_revision = report.content_revision()?;
    report.validate_binding(&report.context_revision)?;
    Ok(report)
}
impl Report {
    fn content_revision(&self) -> Result<String, Error> {
        hash(&(
            (
                &self.engine_version,
                &self.context_revision,
                &self.source_revision,
                self.source_length,
                self.source_feature_count,
                self.source_topology,
                &self.annotations_revision,
                &self.assumptions_revision,
                &self.references_revision,
            ),
            (
                &self.request,
                &self.profile,
                &self.bindings,
                &self.boundaries,
                &self.findings,
                &self.literal_parts,
                &self.literal_searches,
                &self.literal_findings,
                self.literal_search_work,
                &self.statement,
            ),
        ))
    }
    /// Documentary integrity only. Source-based replay is required for science.
    pub fn validate_binding(&self, expected: &str) -> Result<(), Error> {
        if self.bindings.len() > 64
            || self.findings.len() > 256
            || self.literal_parts.len() > 64
            || self.literal_searches.len() > 64
            || self.literal_findings.len() > 64
        {
            return Err(Error::TooLarge);
        }
        // Direct owned callers cannot bypass the canonical location-depth guard.
        let mut locations = collagene_model::bounds::ModelBudget::default();
        for finding in &self.literal_findings {
            if finding.coordinates.len() > 16 {
                return Err(Error::TooLarge);
            }
            for location in &finding.coordinates {
                locations.location(location)?;
            }
        }
        for search in &self.literal_searches {
            if search.matches.len() > 512 {
                return Err(Error::TooLarge);
            }
            for found in &search.matches {
                if found.coordinates.len() > 2 {
                    return Err(Error::TooLarge);
                }
                for location in &found.coordinates {
                    locations.location(location)?;
                }
            }
        }
        bounded(self, MAX_REPORT_BYTES)?;
        self.profile.validate()?;
        self.request.shape()?;
        if self.engine_version != ENGINE_VERSION
            || self.statement != STATEMENT
            || self.profile.reference != self.request.profile
            || self.source_revision != self.request.source_revision
            || self.source_length == 0
            || self.source_length > MAX_BASES
            || self.source_feature_count > 256
            || self.annotations_revision.len() != 64
            || !self
                .annotations_revision
                .bytes()
                .all(|byte| byte.is_ascii_hexdigit() && !byte.is_ascii_uppercase())
            || self.context_revision != expected
            || self.context_revision
                != hash(&(
                    ENGINE_VERSION,
                    super::super::ENGINE_VERSION,
                    crate::sequence::search::SEARCH_ENGINE_VERSION,
                    &self.source_revision,
                    &self.request,
                ))?
            || self.assumptions_revision != hash(&self.request.annotations)?
            || self.references_revision != hash(&self.request.references)?
            || self.report_revision != self.content_revision()?
            || self.bindings.len() != self.request.bindings.len()
            || self.literal_parts.len() != self.request.literal_parts.len()
            || self.literal_search_work > self.request.literal_part_limits.max_work
        {
            return Err(bad("report/source/profile/request integrity changed"));
        }
        let valid_span = |span: &Span| span.start < span.end && span.end <= self.source_length;
        for (node, binding) in self.bindings.iter().zip(&self.request.bindings) {
            if node.binding != *binding
                || node.binding.feature_index >= self.source_feature_count
                || node.interval.as_ref().is_some_and(|span| !valid_span(span))
                || node.interval.is_some() && node.premise != Verdict::Pass
                || node.premise == Verdict::Pass && node.reference_verdict != Verdict::Pass
                || node.reason.len() > 2048
            {
                return Err(bad("invalid bound role evidence"));
            }
        }
        if self.findings.iter().any(|finding| {
            finding.coordinates.len() > 4
                || finding.coordinates.iter().any(|span| !valid_span(span))
                || finding.message.len() > 2048
                || finding.binding_ids.len() > 4
                || finding.feature_indices.len() > 4
                || finding
                    .feature_indices
                    .iter()
                    .any(|index| *index >= self.source_feature_count)
                || finding.source_ids.len() > 8
                || finding.source_ids.iter().any(|id| {
                    !self
                        .profile
                        .definition
                        .sources
                        .iter()
                        .any(|source| &source.id == id)
                })
        }) {
            return Err(bad("invalid bounded vector-class finding"));
        }
        if let Some(boundary) = &self.boundaries {
            if self.request.boundaries.as_ref() != Some(&boundary.selection)
                || boundary.source_spans.len() > 2
                || boundary.source_spans.iter().any(|span| !valid_span(span))
                || boundary.outer_bases > self.source_length
                || boundary.interior_bases > boundary.outer_bases
                || boundary.start_bond > self.source_length
                || boundary.end_bond > self.source_length
            {
                return Err(bad("invalid declared boundary measurement"));
            }
        }
        Ok(())
    }
    pub fn validate_for(&self, source: &Molecule, request: &Request) -> Result<(), Error> {
        self.validate_binding(&revision(source, request)?)?;
        if self != &assess(source, request)? {
            return Err(bad("report differs from current source-based assessment"));
        }
        Ok(())
    }
}
